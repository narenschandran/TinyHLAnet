import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--seed", help = "Seed", required = True,
                    type = int)

parser.add_argument("-f", "--force", action = "store_true",
                    help = "Force rerun")

parser.add_argument("-e", "--embdim", help = "Embedding dimension for Pair Potential", required = True, type = str)
parser.add_argument("-n", "--contacts", help = "Contacts restriction", type = str, default = "allpairs")

parser.add_argument("-p", "--pos_conf", help = "Positional effect configuration", required = False, type = str, default = None)
parser.add_argument("-x", "--effects_conf", help = "Binding independent effects configuration", required = False, type = str, default = None)

parser.add_argument("-d", "--output_dir", help = "Output directory", required = False, type = str, default = None)

parser.add_argument("-t", "--data_type",
                    choices = ["regressand", "binder", "both"],
                    default = "both")

script_path = os.path.realpath(__file__)
script_dir  = os.path.dirname(script_path)
proj_dir    = os.path.join(script_dir, "..")
os.sys.path.append(proj_dir)

from deephlaffylib.Utils import layer_get, layernames
from deephlaffylib.TrainUtils import data_get, model_train

args = parser.parse_args()

if args.effects_conf == "None":
    args.effects_conf = None
if args.pos_conf == "None":
    args.pos_conf = None

from deephlaffylib.DeepHLAffy import DeepHLAffy as Mod
params = {
    'hlalen'      : 44,
    'peplen'      : 9,
    'pp_embdim'   : args.embdim,
    'contacts'    : args.contacts,
    'pos_conf'    : args.pos_conf,
    'effects_conf': args.effects_conf
}

print("-------------------")
print(f"Seed: [{args.seed}]")
for k, v in params.items():
    print(f"{k}: [{v}]")
print("-------------------")

dat = data_get("index", args.data_type)
if args.output_dir is None:
    odir = os.path.join(proj_dir, "models", f"{args.data_type}")
else:
    odir = os.path.join(f"{args.output_dir}", f"{args.data_type}")

hist, perf, mod = model_train(
    Mod, params, dat, args.seed, epochs = 500,
    checkpoints = True, tensorboard = True, base_outdir = odir,
    force = args.force, early_stopping = True)
