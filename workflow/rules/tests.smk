rule make_pipenv:
    output:
        "venv.tar.gz"
    shell:
        tar.using(outputs = ["{output}"])(
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
        "gitignore.done"
    shell:
        tar.using(inputs=["{input}"])(
            "cat {input}/.gitignore && touch {output}"
        )

rule read_pipenv:
    input:
        pip=rules.make_pipenv.output,
        out=rules.read_gitignore.output
    shell:
        boost(
            tar.using(inputs=["{input.pip}"]),
            "cat {input.pip}/pyvenv.cfg"
        )

rule pipenv:
    input: rules.add_gitignore.output

rule test_pipenv_creation:
    shell:
        test_env.script("black -h")

test_script = Pyscript(config["snakemake_dir"], test_env)
rule test_pyscript:
    output:
        **test_script.output(
            foo="out",
        )
    params:
        first="second",
        third=workflow.basedir
    threads: 2
    resources:
        mem_mb=4
    shell:
        test_script(
            "workflow/scripts/test.py",
            params=["third"],
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

