import tensorflow as tf
import keras
from deephlaffylib.Layers import (
        AminoEmbLayer, ANN, ContactsLayer, PosImpBlock,
        EmitOnesMat, EmitZeros,  EmitOnes, EmitSame,

        default_output_layers,

        Int32Cast, DotProduct, MatrixSum, ConcatenateD1, SubsetpHLA
)
from deephlaffylib.Utils import make_uuid

class DeepHLAffy(tf.keras.models.Model):
    '''DeepHLAffy model container

    Args:
        hlalen       : Number of residues in the HLA-I pocket
        peplen       : Number of residues in the peptide
        pp_embdim    : Embedding dimension for the pair potential module
        contacts     : Specific type of contacts assumed by the model
        pos_conf     : Configuration for positional importance module
        effects_conf : Configuration for binding-independent effects module

    Returns:
        Keras model that takes two inputs (HLA and Peptide index+radii for
        individual amino acids) and provides two outputs (binding strength and
        presentation probability).

    Details:
        The "contacts" parameter is configured to take one of the following arguments:
            * allpairs - The model assumes that all contacts are feasible, and
                         does not indicate if certain contacts are more likely.
            * simple   - The model assumes that all contacts observed in
                         crystal structures are feasible and does not
                         indicate if certain contacts are more likely.
            * expdecay - The model assumes that all contacts observed in
                         crystal structures are feasible, but that the feasibilty
                         of contact at a specific residue pair drops
                         exponentially as a function of the specific Cα backbone
                         distance of that residue pair (derived from crystal
                         structures) combined with the Cα centred radii of the
                         specific residues that occupy that location for the
                         given peptide and HLA-I molecule.

        The simplest model, which we term as the baseline model computes
        the binding strength between peptide and HLA-I based purely on the nature
        of the amino acids in the peptide and HLA-I pocket. It assumes that
        all contacts are feasible (allpairs).

        The positional importance module rescales the binding stregnth of any
        given residue pair based on the other residues in the peptide and HLA-I
        pocket. This module is configured using the "pos_conf" parameter. The
        format of the "pos_conf" parameter is a string with the following
        format: <embdim>-<annconf>. The <embdim> should be an integer value
        for configuring the embedding layer for the positional importance. The
        <annconf> should be a standard ANN configuration string (check the
        "ANN" and "parse_dense_conf" functions specified in the
        deephlaffylib.Layers module for more information. Note that the final
        layer in this ANN is supposed to have length equal to hlalen + peplen.
        So in the case where hlalen is 44 and peplen is 9, the ANN configuration
        must end with 53 (44 + 9) nodes.

        In addition to pHLA-I binding strength, binding-independent effects are
        accounted for when calculating presentation probability (qualitative
        output that is fit to immunpeptidome data). The specific Peptide
        involved is assumed to have an additive effect on presentation
        probability. This is meant to represent peptide processing (proteasome
        and TAP transport). This module is configured using the "effects_conf"
        parameter. The configuration of the "effects_conf" is similar to that
        of "pos_conf", but with no restriction on the size of the outputs.

        The "pos_conf" and "effects_conf" parameters accept <None> arguments
        if the model is to be built without these modules. In the case where
        the model has no binding-independent effects, the presentation
        probability is computed only based on the binding strength.
    '''
    def __init__(self, hlalen, peplen, pp_embdim, contacts = "simple",
                 pos_conf = None, effects_conf = None):
        super().__init__()
        self.hlalen      = hlalen
        self.peplen      = peplen
        self.seqlen      = hlalen + peplen
        self.pp_embdim   = int(pp_embdim)

        self.contacts    = contacts
        if self.contacts is None:
            self.contacts = "allpairs"

        self.pos_conf    = pos_conf
        if self.pos_conf == "None":
            self.pos_conf = None

        self.effects_conf = effects_conf
        if self.effects_conf == "None":
            self.effects_conf = None

        # Layers used in the model

        # [1] : Pair potential layers
        self.hla_embmaker       = AminoEmbLayer(self.pp_embdim, name = "hla")
        self.pep_embmaker       = AminoEmbLayer(self.pp_embdim, name = "pep")

        # [2] : Contacts layers
        print(f"--- MODEL INITIALIZATION LOG: Model initialized with contacts of type: [{self.contacts}].")
        self.contactsl    = ContactsLayer(self.contacts, 'contacts', self.seqlen)

        # [3] : Positional importance layers
        if self.pos_conf is None:
            print("--- MODEL INITIALIZATION LOG: No pos_conf given. Model will be built without the positional importance module.")
            self.posimpl = EmitOnesMat(self.seqlen, 'posimp')
        else:
            if '-' not in self.pos_conf:
                raise ValueError("Invalid PosConf")
            print(f"--- MODEL INITIALIZATION LOG: Model initialized with pos_conf: [{self.pos_conf}].")
            pos_embdim, pos_annconf = self.pos_conf.split('-')
            self.posimpl = PosImpBlock(pos_embdim, pos_annconf, name = 'posimp')

        # [4] : Effects layer
        if self.effects_conf is None:
            print("--- MODEL INITIALIZATION LOG: No effects_conf given. Model will be built without the modules for binding-independent effects.")
            #self.pepeffectl = EmitZeros('pepeffect')
            self.pepeffectl = EmitZeros('pepeffect')
        else:
            if '-' not in self.effects_conf:
                raise ValueError("Invalid EffectsConf")
            print(f"--- MODEL INITIALIZATION LOG: Model initialized with effects_conf: [{self.effects_conf}].")
            effects_embdim, effects_annconf = self.effects_conf.split('-')
            effects_embdim = int(effects_embdim)
            self.pepeffectl = keras.Sequential([
                AminoEmbLayer(effects_embdim, name = "pep_effects"),
                ANN(effects_annconf, name = "pep_effects_ann", flatten = True)
            ], name = 'pepeffect')



        # [5] : Output layer
        self.regl, self.bindl = default_output_layers()

        # [6] : Utility layers
        self.dotp        = DotProduct(name = 'dotp')
        self.msum        = MatrixSum(name = 'matsum')
        self.int32_1     = Int32Cast(name = 'i32_1')
        self.int32_2     = Int32Cast(name = 'i32_2')
        self.concat      = ConcatenateD1(name = 'concat')
        self.subset_phla = SubsetpHLA(self.hlalen, name = 'subset_phla')

    def call(self, inputs):
        hlainpfull, pepinpfull = inputs

        # Separate inputs into indices for embeddings and radii for contacts
        hlainp, pepinp = self.int32_1(hlainpfull[:,:,0]), self.int32_2(pepinpfull[:,:,0])
        hlarad, peprad = hlainpfull[:,:,1], pepinpfull[:,:,1]

        # Pair potential - pHLA-I
        hlaemb, pepemb = self.hla_embmaker(hlainp), self.pep_embmaker(pepinp)
        phla_ppmat     = self.dotp(hlaemb, pepemb)

        # Contacts
        contacts      = self.contactsl([hlarad, peprad])
        phla_contacts = self.subset_phla(contacts)

        # Positional importance
        posimp = self.posimpl([hlainp, pepinp])
        phla_posimp = self.subset_phla(posimp)

        # Interaction matrix and binding strength
        phla = phla_ppmat * phla_contacts * phla_posimp
        phla_sum  = self.msum(phla)

        # Binding-independent Effects
        #tapeffect = self.tapl(tapinp)
        pepeffect = self.pepeffectl(pepinp)

        # Model outputs
        # [1] Regression
        reg_intr  = phla_sum
        regout    = self.regl(reg_intr)

        # [2] Binding
        bind_intr = self.concat([phla_sum, pepeffect])
        bindout   = self.bindl(bind_intr)

        return regout, bindout

    @property
    def model_name(self):
        '''Assign unique name to model based on hyperparameter choices

        All models will have the base "deephlaffy", followed by components
        to indicate whether the positional effects module is present ("posmodel"),
        and whether the the binding-independent effects are present ("effects").
        At the end is a unique ID generated based on all the hyperparameters
        of the model instance. This ensures that the final model name is small
        enough to provide module information to humans, but distinct enough
        that models trained with different hyperparameters are not merged
        together.

        The model names generated by this method will be used for segregating
        outputs and performance reports by other functions.
        '''
        base = 'deephlaffy'
        if (self.contacts == 'allpairs') and (self.pos_conf is None) and (self.effects_conf is None):
            base = f"{base}-basemodel"
        else:
            base = f"{base}-{self.contacts}"

        if self.pos_conf is not None:
            base = f"{base}-posmodel"

        if self.effects_conf is not None:
            base = f"{base}-effects"
        
        conf      = self.get_config()
        keys      = sorted(list(conf.keys()))
        comps     = [str(conf[key]) for key in keys]
        comp_str  = "-".join(comps)
        comp_hash = make_uuid(comp_str)
        return f"{base}-{comp_hash}"

    def get_config(self):
        conf = {
            "hlalen": self.hlalen,
            "peplen": self.peplen,
            "pp_embdim": self.pp_embdim,
            "contacts": self.contacts,
            "pos_conf": self.pos_conf,
            "effects_conf": self.effects_conf
        }

        return conf
