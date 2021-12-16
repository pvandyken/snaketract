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
        tar.using(inputs=["{input.pip}"])(
            "cat {input.pip}/pyvenv.cfg"
        )

rule pipenv:
    input: rules.add_gitignore.output

rule test_pipenv_creation:
    shell:
        test_env.script("black -h")

test_script = Pyscript(test_env)
rule test_pyscript:
    output:
        **test_script.output(
            foo="out",
        )
    params:
        first="second",
        third="fourth"
    threads: 2
    resources:
        mem_mb=4
    shell:
        test_script(
            "workflow/scripts/test.py",
            params=["third"],
        )