import toolz as tz

def flatten_args(args: dict) -> list:
    keyword_args = tz.concat(tz.dissoc(args, 'positional_args').items())
    positional_args = args['positional_args']
    return tz.concat([keyword_args, positional_args])