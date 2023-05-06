test_boost = Boost("/tmp", logger)

rule make_pipenv:
    output:
        "venv.tar.gz"
    shell:
        test_boost(
            tar.using(outputs = ["{output}"]),
            "virtualenv {output} && cat {output}/pyvenv.cfg",
        )

rule add_gitignore:
    input:
        rules.make_pipenv.output
    output:
        "gitignore.added"
    shell:
        tar.using(modify=["{input}"])(
            "echo .pyenv >> {input}/.gitignore && touch {output}"
        )

rule read_gitignore:
    input:
        rules.make_pipenv.output
    output:
        touch("gitignore.{i}.done")
    shell:
        test_boost(
            tar.using(inputs=["{input}"]),
            "cat {input}/.gitignore"
        )

rule read_pipenv:
    input:
        pip=rules.make_pipenv.output,
        out=expand(rules.read_gitignore.output, i=range(6))
    shell:
        test_boost(
            tar.using(inputs=["{input.pip}"]),
            "cat {input.pip}/pyvenv.cfg"
        )

rule pipenv:
    input: rules.add_gitignore.output

# rule make_big_tar:
#     input: "data"
#     output: "bigtar.tar.gz"
#     shell:
#         test_boost(
#             tar.using(outputs=["{output}"]),
#             "cp -r data/* {output}"
#         )

# rule read_big_tar:
#     input: "bigtar.tar.gz"
#     output: "read.{i}.txt"
#     shell:
#         test_boost(
#             tar.using(inputs=["{input}"]),
#             "ls {input} |  wc -l > {output}"
#         )

# rule add_to_tar:
#     input: rules.make_big_tar.output
#     output: touch("tar.finished")
#     shell:
#         test_boost(
#             tar.using(modify=["input"]),
#             "touch {input}/foo"
#         )

# rule read_all_big_tars:
#     input:
#         rules.add_to_tar.output,
#         expand(rules.read_big_tar.output, i=range(3)),

rule test_pipenv_creation:
    shell:
        test_env.script("black -h")

test_script = Pyscript(config["snakemake_dir"])
rule test_pyscript:
    output:
        foo="out"
    params:
        first="second",
        third=workflow.basedir
    threads: 2
    resources:
        mem_mb=4
    shell:
        test_env.script(
            test_script(
                "workflow/scripts/test.py",
                params=["third"],
                output=["foo"],
            )
        )

rule write_preamble:
    output: "preamble.{x}.done"
    resources:
        runtime=30,
        mem_mb=1000
    group: "mem_test"
    shell: "touch {output}"

rule consume_memory:
    input: expand("preamble.{x}.done", x=range(5))
    output: "memory.done"
    resources:
        runtime=1,
        mem_mb=1000
    group: "mem_test"
    shell: 
        "touch {output} && (</dev/zero head -c 6G | tail)"
    

rule memory_test:
    input: rules.consume_memory.output

