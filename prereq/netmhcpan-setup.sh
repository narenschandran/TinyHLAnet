SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."

bin_dir="${PROJROOT}/bin"
prereq_dir="${PROJROOT}/prereq"
netmhc_softdir="${prereq_dir}/netmhcpan"

echo "$netmhc_softdir"
f="${netmhc_softdir}/netMHCpan-4.1b.Linux.tar.gz"

if [ ! -f "$f" ]; then
    echo "Unable to detect the netmhcpan archive. Please"
    echo "download it from their website and place it at:"
    echo "[$f]"
    echo "Exiting script..."
    kill -9 $$
fi

data_f="${netmhc_softdir}/data.tar.gz"
if [ ! -f "$data_f" ]; then
    echo "Unable to detect the netmhcpan data. Please"
    echo "download it from their website and place it at:"
    echo "[$data_f]"
    echo "Exiting script..."
    kill -9 $$
fi

# Extract the software
tar -C "$bin_dir" -xvzf "$f"

# Ensure that the main executable is present
netmhc_dir="${bin_dir}/netMHCpan-4.1"
netmhc_f="${netmhc_dir}/netMHCpan"
if [ ! -f "$netmhc_f" ]; then
    echo "Something has gone wrong in extracting the"
    echo "netmhcpan files to:"
    echo "$netmhc_dir"
    echo "Exiting script..."
    kill -9 $$
fi

# Extract data
tar -C "$netmhc_dir" -xvzf "$data_f"

# Customize path to fit this installation
fpath=$(realpath "$netmhc_dir")
tmp_fpath="${fpath}/tmp"
mkdir -p "$tmp_fpath"
sed -i "s|/net/sund-nas.win.dtu.dk/storage/services/www/packages/netMHCpan/4.1/netMHCpan-4.1|$fpath|" "$netmhc_f"
sed -i "s|setenv  TMPDIR  /tmp|setenv  TMPDIR  $tmp_fpath|" "$netmhc_f"
