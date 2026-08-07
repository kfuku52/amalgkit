import os
import re
import shutil
import unicodedata


def cleanup_tmp_amalgkit_files(work_dir='.'):
    try:
        with os.scandir(work_dir) as entries:
            for entry in entries:
                if not entry.name.startswith('tmp.amalgkit.'):
                    continue
                path = entry.path
                try:
                    if entry.is_dir(follow_symlinks=False):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                except FileNotFoundError:
                    continue
    except FileNotFoundError:
        return


_CONTROL_CHARACTER_PATTERN = re.compile(r'[\x00-\x1f\x7f]')
_SPECIES_TOKEN_PATTERN = re.compile(r'^[0-9A-Za-z._-]+$')


def validate_safe_path_component(value, label='Path component'):
    component = str(value).strip()
    if component == '':
        raise ValueError('{} must not be empty.'.format(label))
    if component in ('.', '..'):
        raise ValueError('{} must be a single safe path component: {}'.format(label, component))
    if os.path.isabs(component) or (os.path.basename(component) != component):
        raise ValueError('{} must be a single safe path component: {}'.format(label, component))
    if ('/' in component) or ('\\' in component):
        raise ValueError('{} must not contain path separators: {}'.format(label, component))
    if _CONTROL_CHARACTER_PATTERN.search(component) is not None:
        raise ValueError('{} must not contain control characters.'.format(label))
    return component


def normalize_species_token(scientific_name, label='scientific_name-derived token'):
    raw_name = str(scientific_name).strip()
    if raw_name == '':
        raise ValueError('{} must not be empty.'.format(label))
    if ('/' in raw_name) or ('\\' in raw_name):
        raise ValueError('{} must not contain path separators: {}'.format(label, raw_name))
    if _CONTROL_CHARACTER_PATTERN.search(raw_name) is not None:
        raise ValueError('{} must not contain control characters.'.format(label))
    token = re.sub(r'[^0-9A-Za-z._-]+', '_', raw_name)
    token = re.sub(r'_+', '_', token).strip('_.')
    return validate_safe_path_component(token, label=label)


def resolve_species_token(scientific_name, explicit_token=None, label='species_token'):
    if explicit_token is not None and str(explicit_token).strip() != '':
        token = validate_safe_path_component(explicit_token, label=label)
        if _SPECIES_TOKEN_PATTERN.fullmatch(token) is None:
            raise ValueError(
                '{} must contain only ASCII letters, digits, ".", "_", or "-": {}'.format(
                    label,
                    token,
                )
            )
        return token
    return normalize_species_token(scientific_name, label='scientific_name-derived {}'.format(label))


def validate_unique_species_tokens(species_and_tokens, context='species outputs'):
    species_by_collision_key = {}
    for scientific_name, explicit_token in species_and_tokens:
        normalized_name = str(scientific_name).strip()
        if normalized_name == '':
            continue
        token = resolve_species_token(
            normalized_name,
            explicit_token=explicit_token,
            label='species_token',
        )
        collision_key = unicodedata.normalize('NFC', token).casefold()
        previous = species_by_collision_key.get(collision_key)
        if previous is not None and previous[1] != normalized_name:
            raise ValueError(
                'Scientific names resolve to colliding {} tokens "{}" and "{}": {}, {}'.format(
                    context,
                    previous[0],
                    token,
                    previous[1],
                    normalized_name,
                )
            )
        species_by_collision_key[collision_key] = (token, normalized_name)
    return {
        token: scientific_name
        for token, scientific_name in species_by_collision_key.values()
    }


def build_species_token_map(scientific_names, explicit_tokens=None, context='species outputs'):
    names = [str(value).strip() for value in scientific_names]
    if explicit_tokens is None:
        tokens = [''] * len(names)
    else:
        tokens = ['' if value is None else str(value).strip() for value in explicit_tokens]
    if len(tokens) != len(names):
        raise ValueError('scientific_names and explicit_tokens must have the same length.')
    validate_unique_species_tokens(
        list(zip(names, tokens)),
        context=context,
    )
    token_by_species = {}
    for scientific_name, explicit_token in zip(names, tokens):
        if scientific_name == '':
            continue
        token = resolve_species_token(
            scientific_name,
            explicit_token=explicit_token,
            label='species_token',
        )
        previous = token_by_species.get(scientific_name)
        if previous is not None and previous != token:
            raise ValueError(
                'Scientific name has conflicting species_token values: {} ({}, {})'.format(
                    scientific_name,
                    previous,
                    token,
                )
            )
        token_by_species[scientific_name] = token
    return token_by_species


def _assert_path_is_within_root(path, root, label='Path'):
    path_real = os.path.realpath(path)
    root_real = os.path.realpath(root)
    try:
        common = os.path.commonpath([path_real, root_real])
    except ValueError as exc:
        raise ValueError('{} is outside the allowed root: {}'.format(label, path)) from exc
    if common != root_real:
        raise ValueError('{} is outside the allowed root: {}'.format(label, path))
    return path_real


def ensure_path_within_root(root, path, label='Path'):
    return _assert_path_is_within_root(path=path, root=root, label=label)


def safe_join_component(root, component, label='Path component'):
    component = validate_safe_path_component(component, label=label)
    absolute_root = os.path.abspath(root)
    if os.path.lexists(absolute_root) and os.path.islink(absolute_root):
        raise ValueError('Refusing symbolic-link output root: {}'.format(absolute_root))
    candidate = os.path.join(absolute_root, component)
    if os.path.lexists(candidate) and os.path.islink(candidate):
        raise ValueError('Refusing symbolic-link {}: {}'.format(label, candidate))
    return ensure_path_within_root(
        root=absolute_root,
        path=candidate,
        label=label,
    )


def get_getfastq_run_dir(args, sra_id):
    sra_id = validate_safe_path_component(sra_id, label='Run ID')
    amalgkit_out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(amalgkit_out_dir) and (not os.path.isdir(amalgkit_out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(amalgkit_out_dir))
    getfastq_path = os.path.join(amalgkit_out_dir, 'getfastq')
    if os.path.lexists(getfastq_path) and os.path.islink(getfastq_path):
        raise NotADirectoryError(
            'getfastq path exists but is not a regular directory: {}'.format(getfastq_path)
        )
    getfastq_dir = os.path.realpath(getfastq_path)
    if os.path.exists(getfastq_dir) and (not os.path.isdir(getfastq_dir)):
        raise NotADirectoryError('getfastq path exists but is not a directory: {}'.format(getfastq_dir))
    os.makedirs(getfastq_dir, exist_ok=True)
    run_path = os.path.join(getfastq_dir, sra_id)
    if os.path.lexists(run_path) and os.path.islink(run_path):
        raise NotADirectoryError(
            'Run path exists but is not a regular directory: {}'.format(run_path)
        )
    path_run = _assert_path_is_within_root(
        run_path,
        getfastq_dir,
        label='Run directory',
    )
    if os.path.exists(path_run) and (not os.path.isdir(path_run)):
        raise NotADirectoryError('Run path exists but is not a directory: {}'.format(path_run))
    os.makedirs(path_run, exist_ok=True)
    return path_run
