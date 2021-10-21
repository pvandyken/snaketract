def xvfb_run(config):
    if config.get('x11_srv', False):
        return "xvfb_run"
    return ""