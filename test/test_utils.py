from workflow.lib.utils import flatten_args

def test_args_properly_flattened():
    args = {
        "--preprocessed_data": "test/derivatives/hcp_preproc",
        "--another-arg": "hello world",
        "positional_args": [
            "test",
            "test/derivatives/",
            "participant"
        ]
    }

    flat = flatten_args(args)
    assert list(flat) == [
        "--preprocessed_data",
        "test/derivatives/hcp_preproc",
        "--another-arg",
        "hello world",
        "test",
        "test/derivatives/",
        "participant"
    ]