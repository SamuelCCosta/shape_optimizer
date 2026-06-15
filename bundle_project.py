import tarfile
import os
import glob
import argparse
import subprocess

def create_bundle():
    # Define the core files required for execution
    files_to_include = [
        "shape_optimizer.py",
        "annealing.py",
        "orchestrator.py",
        "db_utils.py",
        "watch_optimization.py",
        "configuration.yaml",
        "configuration-jasminum.yaml"
        "requirements.txt",
        "setup_venv.sh"
    ]
    
    bundle_name = "shape_optimizer_dist.tar.gz"
    
    print(f"Creating bundle: {bundle_name}")
    with tarfile.open(bundle_name, "w:gz") as tar:
        def set_execution_bits(tarinfo):
            # Explicitly set execute permissions for shell scripts
            if tarinfo.name.endswith(".sh"):
                tarinfo.mode |= 0o111 
            return tarinfo

        # Add explicit Python files and config
        for filename in files_to_include:
            if os.path.exists(filename):
                tar.add(filename, filter=set_execution_bits)
                print(f"  [+] {filename}")
        
        # Specifically capture Linux shared objects for the solver
        for match in glob.glob("square_solver*.so"):
            tar.add(match, filter=set_execution_bits)
            print(f"  [+] {match}")

    return bundle_name

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bundle the shape optimizer project and optionally send via SSH.")
    parser.add_argument("--ssh", help="Remote destination (e.g., user@host:/path/to/dir/)")
    parser.add_argument("-p", "--port", help="SSH port", type=str)
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output for SCP")
    parser.add_argument("--force-password", action="store_true", help="Force password authentication (ignore local keys)")
    args = parser.parse_args()

    bundle_name = create_bundle()

    if args.ssh:
        print(f"\nTransferring '{bundle_name}' to {args.ssh}...")
        scp_cmd = ["scp"]
        if args.verbose:
            scp_cmd.append("-v")
        if args.port:
            scp_cmd.extend(["-P", args.port])
        if args.force_password:
            # These options tell SSH to skip key-based auth and not use the ssh-agent
            scp_cmd.extend([
                "-o", "PreferredAuthentications=password",
                "-o", "PubkeyAuthentication=no"
            ])
        scp_cmd.extend([bundle_name, args.ssh])

        result = subprocess.run(scp_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print("Transfer successful.")
        else:
            print(f"Transfer failed with error:\n{result.stderr}")
    else:
        print(f"\nBundle complete. You can now transfer '{bundle_name}' to your target environment.")