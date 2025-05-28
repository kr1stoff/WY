def generate_shell(output_shell, content, outshell, finish_string=None):
    if finish_string is None:
        finish_string = "Still_waters_run_deep"

    # Remove existing output_shell files
    existing_files = glob.glob(f"{output_shell}.*")
    for file in existing_files:
        os.remove(file)

    # Write content to the output_shell file
    with open(output_shell, "w") as out_file:
        out_file.write("#!/bin/bash\n")
        out_file.write("echo ==========start at : `date` ==========\n")
        out_file.write("set -e \n")
        out_file.write(f"{content} \n")
        out_file.write("echo ==========end at : `date` ========== && \n")
        out_file.write(f"echo {finish_string} 1>&2 && \\\n")
        out_file.write(f"echo {finish_string} > {output_shell}.sign\n")

    run_command = f"bash {output_shell} 1>{output_shell}.e 2>{output_shell}.o"
    outshell.append(run_command)
