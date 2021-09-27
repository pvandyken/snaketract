import itertools as it
import toolz as tz

def flatten_args(args: dict) -> list:
    keyword_args = it.chain(tz.dissoc(args, 'positional_args').items())
    positional_args = args['positional_args']
    return it.chain([keyword_args, positional_args])