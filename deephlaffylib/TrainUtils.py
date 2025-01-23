import os
import re
import time

import numpy as np
import pandas as pd

import tensorflow as tf
import keras
from keras.callbacks import Callback, EarlyStopping, ModelCheckpoint, TensorBoard

from nbslpy._utils import savePickle, readPickle
from deephlaffylib.Losses import custom_mse, custom_bce
from deephlaffylib.Utils import (
        model_load, writeDict, readDict,
        MSE, BCE, SRCC, AUC, PRAUC
)

#------------------------------------------------------------------------------#
#              Find and load preprocessed data for model training              #
#------------------------------------------------------------------------------#

# This script is hardcoded to expect preprocessed datasets in the "datasets"
# directory accessed from the project root directory.
script_dir    = os.path.dirname(__file__)
projroot      = os.path.join(script_dir, "..")
base_data_dir = os.path.join(projroot, "datasets")


def pklfiles_get(d):
    '''Retrieve the processed pkl files and split them by split and preprocessing type.'''
    fs = [f for f in os.listdir(d) if f[-3:] == "pkl"]

    splits = set([f.split(".")[0] for f in fs])
    types = set([f.split(".")[1] for f in fs])

    fs_dict = {}
    for tp in types:
        fs_dict[tp] = {}
        for sp in splits:
            fs_dict[tp][sp] = os.path.join(d, f"{sp}.{tp}.pkl")

    return fs_dict


def data_assemble(pkl_file, inpo, outo):
    '''Assembles input and output datapoints into a two element list.

    The inputs are placed first in the list followed by the outputs. The order
    of the features in both elements are specified through "inpo" and "outo"
    respectively.
    '''
    def _data_assemble(d, o):
        return [d[oo] for oo in o]

    dat = readPickle(pkl_file)
    return [_data_assemble(dat, inpo), _data_assemble(dat, outo)]


def data_get(data_type, subset="both"):
    '''Wrapper function to retrieve data in a format ready to be fed for model training.'''
    inpo = ["hla", "pep"]
    outo = ["regressand", "binder"]
    data = {}

    # The function is hardcoded to expect the preprocessed datasets to
    # be present in the following location:
    proc_dir = os.path.join(base_data_dir, "proc")

    # We retrieve the preprocessed datasets from the following location
    # and ensure that only the required data (regerssion/binder/both) is
    # produced to the user.
    data_files = pklfiles_get(proc_dir)
    fs = data_files[data_type]

    for sp in fs:
        dat = data_assemble(fs[sp], inpo, outo)
        if subset != "both":
            if subset == "regressand":
                l = dat[1][0] > (-1.0)
            elif subset == "binder":
                l = dat[1][1] > (-1.0)
            else:
                raise ValueError(
                    "Data subset should be 'both', 'regressand' or 'binder'"
                )
            ind = np.argwhere(l)[:, 0]
            dat[0] = [inp[ind, ...] for inp in dat[0]]
            dat[1] = [inp[ind, ...] for inp in dat[1]]
        data[sp] = dat

    return data

#------------------------------------------------------------------------------#
#                  Performance assessement for dataset splits                  #
#------------------------------------------------------------------------------#
def dataset_perf(mod, dataset):
    """
    Wrapper function to compute different error statistics (both for
    regression and label outputs) for a given model.
    """
    inp, out = dataset
    pred = mod.predict(inp)
    perf = {
        "SRCC": SRCC(out[0], pred[0]),
        "MSE": MSE(out[0], pred[0]),
        "BCE": BCE(out[1], pred[1]),
        "AUC": AUC(out[1], pred[1]),
        "PRAUC": PRAUC(out[1], pred[1]),
    }
    return perf


def perf_report_generate(mod, data, seed):
    """
    Generates performance reports for training-val-test splits and provides
    the worst-case scenario for all cases.
    """
    perfs = {}
    for split in ["train", "val", "test"]:
        perfs[split] = dataset_perf(mod, data[split])

    perf_meas = list(perfs["train"].keys())

    split_ord = ["train", "val", "test"]
    vlist = []

    # Decides whether 'worst-case' means minimum or maximum value for each
    # of the summary statistics.
    for meas in perf_meas:
        vals = []
        for split in split_ord:
            vals.append(perfs[split][meas])
        if meas in ["SRCC", "AUC", "PRAUC"]:
            vals.append(min(vals))
        elif meas in ["MSE", "BCE"]:
            vals.append(max(vals))
        else:
            raise ValueError(f"Unknown measure: [{meas}]")
        vlist.append(np.asarray([meas] + vals))

    perf_df = pd.DataFrame(vlist)
    perf_df.columns = ["measure", "train", "val", "test", "worst-case"]
    perf_df.loc[:, "seed"] = seed
    perf_df.loc[:, "model_id"] = mod.model_name

    return perf_df


#------------------------------------------------------------------------------#
#                               Custom Callbacks                               #
#------------------------------------------------------------------------------#
def DefaultEarlyStopping():
    """Early stopping on validation loss to prevent overfitting"""
    return EarlyStopping(
        min_delta=0.001,
        patience=15,
        monitor="val_loss",
        start_from_epoch=15,
        restore_best_weights=True,
    )

class PerfReportAtEachEpoch(Callback):
    """Generate a perf-report on input data at the end of every epoch"""
    def __init__(self, mod, data, seed, odir):
        super().__init__()
        self.mod = mod
        self.data = data
        self.seed = seed
        self.odir = odir
        if not os.path.exists(self.odir):
            os.makedirs(self.odir, exist_ok = True)

    def on_epoch_end(self, epoch, logs = None):
        perf = perf_report_generate(self.mod, self.data, self.seed)
        perf_file = os.path.join(self.odir, f"perf-e{epoch:04d}.tsv")
        perf.to_csv(perf_file, sep = '\t', index = False)


#------------------------------------------------------------------------------#
#                                Model training                                #
#------------------------------------------------------------------------------#
# We have a reasonably sophisticated system to ensure that model training is   #
# not repeat if we can detected instances of prior training. Our system also   #
# provides the options to track checkpoints and use Tensorboard using the      #
# callbacks provided by Tensorflow. We also configure a default early stopping #
# callback (based on the callback provided by Tensorflow) to prevent           #
# overfitting.                                                                 #
#                                                                              #
# Out intention with this system is to use EarlyStopping to figure out the     #
# best hyperparameter set in a time-efficient manner, followed by running      #
# the best hyperparameter set without early stopping for a large number of     #
# epochs, where we track the performance across all epochs.                    #
#------------------------------------------------------------------------------#

def best_model_from_chkpts(chkpt_dir):
    '''Gets the best model from a checkpoints directory based on the loss value.'''
    def parse_fname(f):
        bname = os.path.basename(f)
        fbase = re.sub("[.]keras$", "", bname)
        epoch, loss = fbase.split("-")
        return epoch, float(loss)

    fs = [os.path.join(chkpt_dir, fname) for fname in os.listdir(chkpt_dir)]
    fs.sort()

    epochs = []
    losses = []
    fpaths = []

    for f in fs:
        epoch, loss = parse_fname(f)
        epochs.append(epoch)
        losses.append(loss)

    datf = pd.DataFrame({"epochs": epochs, "loss": losses, "f": fs})
    return datf.sort_values("loss").iloc[0].loc["f"]



def model_train(modfn, params, data, seed, base_outdir=None,
                epochs=500, verbose=2, early_stopping=False,
                checkpoints=False, tensorboard=False, force=False,
                mse_scale=1.0, bce_scale=1.0, report_all = False):
    '''
    Wrapper function to set up model training

    Args:
        modfn          : The model class (used to instantiate the model).
        params         : Arguments passed to the model class.
        data           : Dictionary with train-val-test splits of the dataset.
        base_outdir    : Base output directory.
        epochs         : Maximum number of epochs used for the training process.
        verbose        : Control verbosity during model fitting.
        early_stopping : Enable/disable early stopping to prevent overfitting.
        checkpoints    : Enable/disable saving model checkpoints. If enabled,
                         the model that is reported in the output directory
                         will be the one with the best performance.
        tensorboard    : Enable/disable production of Tensorboard reports
        force          : Perform model training from scratch even if prior
                         model exists.
        mse_scale      : Scaling factor for regression loss function.
        bce_scale      : Scaling factor for binder label loss function.
        report_all     : Compute performance report at the end of each epoch.
    '''
    def data_check(data):
        splits = ["train", "val", "test"]
        for split in splits:
            if split not in data.keys():
                raise ValueError(f"The [{split}] split is not present")

        return data["train"], data["val"], data["test"]

    # Check [1]: Ensure that the train, val and test splits
    # are present in the input data.
    data_check(data)
    traininp, trainout = data["train"]
    val = data["val"]

    if base_outdir is not None:
        mdl_name = modfn(**params).model_name
        outdir = os.path.join(base_outdir, mdl_name, f"seed-{seed:02d}")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        perf_file  = os.path.join(outdir, "perf.tsv")
        conf_file  = os.path.join(outdir, "conf.dct")
        model_file = os.path.join(outdir, "model.keras")
        hist_file  = os.path.join(outdir, "hist.pkl")
        time_file  = os.path.join(outdir, "time.txt")

        if (not force) and (os.path.exists(perf_file) and os.path.exists(model_file) and os.path.exists(conf_file) and os.path.exists(hist_file)):
            mod = model_load(modfn, model_file)
            hist = readPickle(hist_file)
            perf = pd.read_csv(perf_file, sep = '\t')

            print("")
            print("-----------------------------------------------------------")
            print("The existence of the following files indicates the presence")
            print("of a pretrained model. If you want to do a fresh training")
            print("delete these files and continue:")
            print(perf_file)
            print(model_file)
            print(conf_file)
            print("-----------------------------------------------------------")
            print("")

            pcols = ['measure', 'train', 'val', 'test', 'worst-case', 'seed']
            print(perf.loc[:,pcols])
            return hist, perf, mod
    else:
        outdir = None

    callbacks = []
    if early_stopping:
        callbacks.append(DefaultEarlyStopping())

    if checkpoints and (base_outdir is None):
        raise ValueError("Specify base_outdir if you want to store checkpoints")

    if tensorboard and (base_outdir is None):
        raise ValueError("Specify base_outdir if you want to use Tensorboard")

    if checkpoints and (base_outdir is not None):
        chkptdir = os.path.join(outdir, "chkpt")
        if not os.path.exists(chkptdir):
            os.makedirs(chkptdir)
        chkptfile = os.path.join(chkptdir, "{epoch:04d}-{val_loss:0.4f}.keras")
        callbacks.append(ModelCheckpoint(chkptfile))

    if tensorboard and (base_outdir is not None):
        logdir = os.path.join(outdir, "logs")
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        callbacks.append(TensorBoard(logdir, histogram_freq=1))

    # Set seed once before model initialization
    tf.keras.utils.set_random_seed(seed)
    mod = modfn(**params)
    mod.compile(
        optimizer="adam", loss=[custom_mse(mse_scale), custom_bce(bce_scale)]
    )

    if report_all:
        report_odir = os.path.join(outdir, 'perf-by-epoch')
        callbacks.append(PerfReportAtEachEpoch(mod, data, seed, report_odir))

    # Setting seed once again before training starts for
    # the sake of safety.
    start = time.time()
    tf.keras.utils.set_random_seed(seed)
    hist = mod.fit(
        traininp,
        trainout,
        validation_data=val,
        epochs=epochs,
        callbacks=callbacks,
        verbose=verbose,
    )
    end = time.time()
    elaps = end - start
    time_str = f"Model training time: {elaps:0.2f} seconds"
    print(time_str)
    if base_outdir is not None:
        with open(time_file, 'w') as f:
            f.writelines(time_str + '\n')

    if checkpoints:
        best_model_f = best_model_from_chkpts(chkptdir)
        print(f"Selecting best model from checkpoints:\n[{best_model_f}]")
        mod = model_load(modfn, best_model_f)

    perf = perf_report_generate(mod, data, seed)
    conf = mod.get_config()

    if outdir is not None:
        mod.save(model_file)
        writeDict(perf, perf_file)
        perf.to_csv(perf_file, sep = '\t', index = False)
        writeDict(conf, conf_file)
        savePickle(hist, hist_file)


        # Just testing that it can load the model
        tmp = model_load(modfn, model_file)
        del tmp
    pcols = ['measure', 'train', 'val', 'test', 'worst-case', 'seed']
    print(perf.loc[:,pcols])
    return hist, perf, mod
