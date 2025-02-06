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
    local data_type='regressand'
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

pos_confs_gen() {
    actvs1='linear tanh gelu relu sigmoid'
    actvs2='relu sigmoid'
    for actv1 in $(echo "$actvs1"); do
    for actv2 in $(echo "$actvs2"); do
        echo "16-${actv1}16_${actv2}53"
    done
    done
}

# These hyperparameters are fixed based on the prior experiments
contact_type='simple'
embdim='16'

# We first optimize the activations
pos_confs=$(pos_confs_gen)
for pos_conf in $(echo "$pos_confs"); do
    run_model "$embdim" "$contact_type" "$pos_conf"
done

# Then we optimize the embedding dimension and hidden nodes in the ANN
dims='256 128 64 32 16 8 4'
for pos_dim in $(echo "$dims"); do
for pos_node in $(echo "$dims"); do
    run_model "$embdim" "$contact_type" "${pos_dim}-sigmoid${pos_node}_sigmoid53"
done
done
