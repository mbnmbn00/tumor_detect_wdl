version 1.0

workflow variant_calling {
	input {
		File read1_file
    File read2_file
    String qced1_file
    String qced2_file
	}
  call run_fastp {
    input:
      read1_file = read1_file,
      read2_file = read2_file,
      qced1_file = qced1_file,
      qced2_file = qced2_file
  }
}

task run_fastp {
	input {
		File read1_file
    File read2_file
    String qced1_file
    String qced2_file
	}
	command {
    fastp \
      --in1 ${read1_file} --out1 ${qced1_file} \
      --in2 ${read2_file} --out2 ${qced2_file}
  }
	output {
		String out = read_string(stdout())
	}
  runtime {
    docker: "tumor_img:latest"
  }
}
