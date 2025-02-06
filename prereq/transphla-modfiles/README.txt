We fixed two small bugs in the pHLAIformer interface:
# Fix indentation error
sed -i 's/        log = Logger(errLogPath)/    log = Logger(errLogPath)/' pHLAIformer.py

# Fix bug in iteration variable (hla->hla_seq)
sed -i 's/    if not (pep.isalpha() and hla.isalpha()):/    if not (pep.isalpha() and hla_seq.isalpha()):/' pHLAIformer.py

# We also changed all file paths so that the code can be run from anywhere. The changed lines are:

model.py:2:SCRIPT_DIR = os.path.dirname(__file__)
model.py:5:sys.path.append(SCRIPT_DIR)
model.py:58:vocab = np.load(os.path.join(SCRIPT_DIR, 'vocab_dict.npy'), allow_pickle = True).item()

mutation.py:2:SCRIPT_DIR = os.path.dirname(__file__)
mutation.py:5:sys.path.append(SCRIPT_DIR)
mutation.py:51:        aatype_position = np.load(os.path.join(SCRIPT_DIR, './Attention/peptideAAtype_peptidePosition/{}_Length{}.npy'.format(hla, length)),
mutation.py:54:        aatype_position_num = np.load(os.path.join(SCRIPT_DIR, './Attention/peptideAAtype_peptidePosition_NUM/{}_Length{}_num.npy'.format(hla, length)),
mutation.py:58:        aatype_position = np.load(os.path.join(SCRIPT_DIR, './Attention/peptideAAtype_peptidePosition/Allsamples_Alllengths.npy'),
mutation.py:60:        aatype_position_num = np.load(os.path.join(SCRIPT_DIR, './Attention/peptideAAtype_peptidePosition_NUM/Allsamples_Alllengths_num.npy'),


pHLAIformer.py:2:SCRIPT_DIR = os.path.dirname(__file__)
pHLAIformer.py:5:sys.path.append(SCRIPT_DIR)
pHLAIformer.py:166:model_file = os.path.join(SCRIPT_DIR, 'pHLAIformer.pkl')
