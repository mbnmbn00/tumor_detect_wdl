version development

# Define the inputs for the workflow
workflow runPythonScripts {
  Array[File] read_files
  File assembly_file
  String host_path

  # Define the call for running the first Python script (run_fastp.py)
  call run_fastp {
    input:
      read_files = read_files,
      assembly_file = assembly_file
      docker_image = tumor_img
      host_path = host_path
  }

  # Define the call for running the second Python script (run_minimap2.py)
  call run_minimap2 {
    input:
      read_files = run_fastp.output.fastq_output,
      assembly_file = assembly_file
      docker_img = tumor_img
      host_path = host_path
  }

  # Define the call for running the third Python script (run_mutect2.py)
  call run_mutect2 {
    input:
      bam_file = run_minimap2.output.bam_output
      docker_img = tumor_img
      host_path = host_path
  }
}

# Define the tasks for running the Python scripts
task run_fastp {
  Array[File] read_files
  File assembly_file
  command {
    python run_fastp.py --read_files ${sep=" " read_files} --assembly_file ${assembly_file}
  }
  output {
    File fastq_output = "fastp_output.fastq"
  }
  runtime {
    docker: "tumor_img",
    volumes: [
      host_path:/data
    ]
  }
}

task run_minimap2 {
  Array[File] read_files
  File assembly_file
  command {
    python run_minimap2.py --read_files ${sep=" " read_files} --assembly_file ${assembly_file}
  }
  output {
    File bam_output = "minimap2_output.bam"
  }
  runtime {
    docker: "tumor_img"
  }
}

task run_mutect2 {
  File bam_file
  command {
    python run_mutect2.py --bam_file ${bam_file}
  }
  output {
    File vcf_output = "mutect2_output.vcf"
  }
  runtime {
    docker: "tumor_img"
  }
}
