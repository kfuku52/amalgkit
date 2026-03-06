from types import SimpleNamespace


def namespace_to_dict(args):
    data = {}
    if hasattr(args, '__dict__'):
        data.update(vars(args))
    for name in dir(args):
        if name.startswith('__'):
            continue
        if name in data:
            continue
        try:
            value = getattr(args, name)
        except AttributeError:
            continue
        if callable(value):
            continue
        data[name] = value
    return data


def clone_namespace(args, **overrides):
    data = namespace_to_dict(args)
    data.update(overrides)
    return SimpleNamespace(**data)
