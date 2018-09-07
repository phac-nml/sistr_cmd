import logging
import numpy as np


def listattrs(x):
    """Get all instance and class attributes for an object
    
    Get all instance and class attributes for an object except those that start
    with "__" (double underscore).
    
    __dict__ of an object only reports the instance attributes while dir() 
    reports all of the attributes of an object including private ones.
    
    Callable attrs are filtered out.
    
    Args:
        x (object): Some object
    
    Returns:
        list str: List of non-callable non-private attributes of object x
    """
    return [attr for attr in dir(x) if not attr.startswith("__") and not callable(getattr(x, attr))]


def to_dict(x, depth, exclude_keys=set(), depth_threshold=8):
    """Transform a nested object/dict/list into a regular dict
    
    json.dump(s) and pickle don't like to un/serialize regular Python objects so
    this function should handle arbitrarily nested objects to be serialized to
    regular string, float, int, bool, None values.
    
    This is a recursive function so by default it will exit at a certain depth (depth_threshold=8).
    
    Args:
        x (object): Some object to dict-ify unless x is a scalar/literal then return x as is
        depth (int): Starting depth must be 0 (cannot supply default value due to weird Pythonisms)
        exclude_keys (set): Keys to avoid adding to the output dict
        depth_threshold (int): object/dict nesting depth to stop at
    
    Returns:
        dict: dict with only scalar/literal leaf values    
    """
    if x is None or isinstance(x, (str, int, float, bool)):
        return x
    if isinstance(x, np.int_):
        return int(x)
    if isinstance(x, np.int64):
        return int(x)
    if isinstance(x, np.float_):
        return float(x)
    if isinstance(x, np.float64):
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if depth + 1 > depth_threshold: return {}
    if isinstance(x, list):
        out = []
        for v in x:
            tmp = to_dict(v, depth + 1, exclude_keys, depth_threshold)
            if tmp == {}: continue
            out.append(tmp)
        return out
    out = {}
    if isinstance(x, dict):
        for k, v in x.items():
            if k in exclude_keys: continue
            if not isinstance(k, (str,)):
                k = str(k)
            tmp = to_dict(v, depth + 1, exclude_keys, depth_threshold)
            if tmp == {}: continue
            out[k] = tmp
        return out
    for attr in listattrs(x):
        if attr in exclude_keys: continue
        v = getattr(x, attr)
        tmp = to_dict(v, depth + 1, exclude_keys, depth_threshold)
        if tmp == {}: continue
        out[attr] = tmp
    return out


def _recur_flatten(key, x, out, sep='.'):
    """Helper function to flatten_dict
    
    Recursively flatten all nested values within a dict
    
    Args:
        key (str): parent key
        x (object): object to flatten or add to out dict
        out (dict): 1D output dict
        sep (str): flattened key separator string
    
    Returns:
        dict: flattened 1D dict
    """
    if x is None or isinstance(x, (str, int, float, bool)):
        out[key] = x
        return out
    if isinstance(x, list):
        for i, v in enumerate(x):
            new_key = '{}{}{}'.format(key, sep, i)
            out = _recur_flatten(new_key, v, out, sep)
    if isinstance(x, dict):
        for k, v in x.items():
            new_key = '{}{}{}'.format(key, sep, k)
            out = _recur_flatten(new_key, v, out, sep)
    return out


def flatten_dict(x):
    """Flatten a dict
    
    Flatten an arbitrarily nested dict as output by to_dict
    
    .. note:: 
    
        Keys in the flattened dict may get very long. 
        
    Args:
        x (dict): Arbitrarily nested dict (maybe resembling a tree) with literal/scalar leaf values
    
    Returns:
        dict: flattened 1D dict    
    """
    out = {}
    for k, v in x.items():
        out = _recur_flatten(k, v, out)
    return out


def write_json(fh, output):
    import json
    json.dump(output, fh)


def write_pickle(fh, output):
    import pickle
    pickle.dump(output, fh)


def write_csv(fh, output):
    import pandas as pd
    df = pd.DataFrame(output)
    df.to_csv(fh, index=False)


def write_tab(fh, output):
    import pandas as pd
    df = pd.DataFrame(output)
    df.to_csv(fh, index=False, sep='\t')


fmt_to_write_func = {'json': write_json,
                     'pickle': write_pickle,
                     'csv': write_csv,
                     'tab': write_tab, }


def write(dest, fmt, serovar_predictions, more_results=0):
    assert isinstance(serovar_predictions, list)
    if not fmt in fmt_to_write_func:
        logging.warning('Invalid output format "%s". Defaulting to "json"', fmt)
        fmt = 'json'
    if '.' + fmt not in dest:
        dest += '.' + fmt
    logging.info('Writing output "%s" file to "%s"', fmt, dest)
    fh = open(dest, 'w')
    try:
        # write in whatever format necessary
        write_func = fmt_to_write_func[fmt]
        exclude_keys_in_output = {'blast_results', 'sseq'}
        if more_results >= 2:
            exclude_keys_in_output.remove('blast_results')
            exclude_keys_in_output.remove('sseq')
        elif more_results == 1:
            exclude_keys_in_output.remove('sseq')
        if fmt in {'pickle', 'json'}:
            output_dict = [to_dict(v, 0, exclude_keys=exclude_keys_in_output) for v in serovar_predictions]
        else:
            if more_results > 0:
                output_dict = [flatten_dict(to_dict(v, 0, exclude_keys=exclude_keys_in_output)) for v in
                               serovar_predictions]
            else:
                output_dict = [to_dict(v, 0, exclude_keys=exclude_keys_in_output, depth_threshold=1) for v in
                               serovar_predictions]
        write_func(fh, output_dict)
    finally:
        fh.close()
