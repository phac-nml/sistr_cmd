import logging
from src.serovar_prediction import SerovarPrediction


def write_json(fh, output):
    import json
    json.dump(output, fh)

def write_pickle(fh, output):
    import cPickle
    cPickle.dump(output, fh)

def write_csv(fh, output):
    import pandas as pd
    df = pd.DataFrame([output])
    df.to_csv(fh, index=False)

fmt_to_write_func = {'json': write_json,
                     'pickle': write_pickle,
                     'csv': write_csv}


def write(dest, fmt, serovar_prediction):
    assert isinstance(serovar_prediction, SerovarPrediction)
    if not fmt in fmt_to_write_func:
        logging.warn('Invalid output format "{}". Defaulting to "json"'.format(fmt))
        fmt = 'json'

    if '.' + fmt not in dest:
        dest += '.' + fmt


    logging.info('Writing output "{}" file to "{}"'.format(fmt, dest))
    fh = open(dest, 'w')
    try:
        # write in whatever format necessary
        write_func = fmt_to_write_func[fmt]
        if fmt == 'pickle':
            output_dict = serovar_prediction.__dict__
        else:
            output_dict = {}
            for k,v in serovar_prediction.__dict__.items():
                if isinstance(v, (str, float, int)):
                    output_dict[k] = v
        write_func(fh, output_dict)
    except Exception as ex:
        logging.error(ex)
    finally:
        fh.close()

