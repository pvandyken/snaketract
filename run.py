#!/usr/bin/env python3
import os
from snakebids.app import SnakeBidsApp
import snakebids



def main():
    app = SnakeBidsApp(os.path.abspath(os.path.dirname(__file__)))
    app.config["output_dir"] = os.path.join(app.config["output_dir"], "prepdwi_recon")
    app.run_snakemake()


if __name__ == "__main__":
    main()