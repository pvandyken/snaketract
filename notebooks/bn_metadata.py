import pandas as pd
import copy

def read_metadata():
    return (
        pd.concat([
            copy.copy(
                orig := (
                    raw := pd.read_csv("resources/brainnetome-regions.csv")
                )
                .assign(
                    **dict(
                        raw["Modified cyto-architectonic"]
                        .str.split(r',\s?', expand=True)
                        .rename({0: "Name", 1: "Long Name"}, axis=1)
                        .drop(columns=[2])
                    ),
                    **dict(
                        raw["Gyrus"]
                        .str.split(r',\s?', expand=True)
                        .rename({0: "Gyrus Abbr", 1: "Gyrus"}, axis=1)
                    )
                )
                .drop(
                    columns=[
                        "Modified cyto-architectonic",
                        "Left and right hemispheres"
                    ]
                )
            )
            .drop(columns=["rh.MNI (X, Y, Z)", "Label ID.R"])
            .rename(columns={"lh.MNI (X,Y,Z)": "MNI", "Label ID.L": "Label ID"})
            .assign(hemisphere="R"),

            orig
            .drop(columns=["lh.MNI (X,Y,Z)", "Label ID.L"])
            .rename(columns={"rh.MNI (X, Y, Z)": "MNI", "Label ID.R": "Label ID"})
            .assign(hemisphere="L")
        ])
        .astype({"Label ID": int})
        .assign(
            Lobe=lambda df: df["Lobe"].str.strip()
        )
        .set_index("Label ID")
        .sort_index()
        # .sort_values(["hemisphere"])
    )