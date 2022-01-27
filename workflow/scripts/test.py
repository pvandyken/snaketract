from snakeboost import snakemake_args
from pathlib import Path

args = snakemake_args()

out = Path("out")

with out.open("w") as f:
    f.write(
        f"""input={args.input}
        output={args.output}
        params={args.params}
        resources={args.resources}
        threads={args.threads}
        log={args.log}
        """
    )
fi = 160_000_00