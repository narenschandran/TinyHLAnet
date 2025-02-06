# Script setup to activate environment and to select the correct
# Python version.
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"
MDL_DIR="${PROJROOT}/models/deephlaffy"

env_dir="${PROJROOT}/.env/deephlaffy"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi
PYF="${SCRIPT_DIR}/train.py"

run_model() {
    # Fixed hyperparameters
    local data_type='both'
    local output_dir="${MDL_DIR}"
    local nmodels=3

    local embdim="$1"
    local contact_type="$2"
    local pos_conf="$3"
    local effects_conf="$4"

    if [ -z "$pos_conf" ]; then
        local pos_conf="None"
    fi

    if [ -z "$effects_conf" ]; then
        local effects_conf="None"
    fi

    seedlim=$(echo "${nmodels} - 1" | bc)
    for seed in $(seq 0 "$seedlim"); do
    echo "$PY" "$PYF" -t "$data_type" \
        -e "$embdim" -n "$contact_type" \
        -p "$pos_conf" -x "$effects_conf" \
        -d "$output_dir" -s "$seed"

    # If not a dry run, execute the model training
    if [ -z "$DRYRUN" ]; then
        "$PY" "$PYF" -t "$data_type" \
            -e "$embdim" -n "$contact_type" \
            -p "$pos_conf" -x "$effects_conf" \
            -d "$output_dir" -s "$seed"
    fi
    done
}

# These hyperparameters are fixed based on prior experiments.
contact_type='simple'
embdim='16'
pos_conf='16-sigmoid128_sigmoid53'

# We fit the models based on just binding affinity to the full dataset
# (quantitative + qualitative data. The use of the full dataset is specified
# within the run_model function in this script).
run_model "$embdim" 'allpairs'
run_model "$embdim" "$contact_type"
run_model "$embdim" "$contact_type" "$pos_conf"

# This is the best pos_conf
# Next, we fit models with binding-independent effects
effects_conf='16-linear16_linear1'
run_model "$embdim" "$contact_type" "$pos_conf" "$effects_conf"

# Optimize activation for binding-independent effects
in_actvs='linear relu sigmoid tanh gelu'
out_actv='linear'
for in_actv in $(echo "$in_actvs"); do
    run_model "$embdim" "$contact_type" "$pos_conf" "16-${in_actv}16_${out_actv}1"
done

# Optimize embedding dimension and hidden node count
dims='256 128 64 32 16 8 4'
nodes='256 128 64 32 16 8 4'

for effects_emb in $(echo "$dims"); do
for effects_node in $(echo "$nodes"); do
    run_model "$embdim" "$contact_type" "$pos_conf" "${effects_emb}-gelu${effects_node}_linear1"
done
done
