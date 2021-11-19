rule make_pipenv:
    output:
        "venv.tar.gz"
    shell:
        tar(
            outputs = ["{output}"],
            cmd="virtualenv {output} && cat {output}/pyvenv.cfg",
        )

rule read_gitignore:
    input:
        rules.make_pipenv.output
    output:
        "gitignore.done"
    shell:
        tar(
            inputs=["{input}"],
            cmd="cat {input}/.gitignore && touch {output}"
        )

rule add_gitignore:
    input:
        rules.make_pipenv.output
    output:
        "gitignore.added"
    shell:
        tar(
            modify=["{input}"],
            cmd="echo .pyenv >> {input}/.gitignore && touch {output}"
        )

rule read_pipenv:
    input:
        pip=rules.make_pipenv.output,
        out=rules.read_gitignore.output
    shell:
        tar(
            inputs=["{input.pip}"],
            cmd="cat {input.pip}/pyvenv.cfg"
        )

rule pipenv:
    input: rules.add_gitignore.output

