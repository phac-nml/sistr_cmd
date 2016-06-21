import logging


def write_json(fh, output):
    import json
    json.dump(output, fh)

def write_pickle(fh, output):
    import cPickle
    cPickle.dump(output, fh)

def write_csv(fh, output):
    import pandas as pd
    df = pd.DataFrame(output)
    df.to_csv(fh, index=False)

fmt_to_write_func = {'json': write_json,
                     'pickle': write_pickle,
                     'csv': write_csv}


def write(dest, fmt, serovar_predictions):
    assert isinstance(serovar_predictions, list)
    if not fmt in fmt_to_write_func:
        logging.warn('Invalid output format "%s". Defaulting to "json"', fmt)
        fmt = 'json'

    if '.' + fmt not in dest:
        dest += '.' + fmt


    logging.info('Writing output "%s" file to "%s"', fmt, dest)
    fh = open(dest, 'w')
    try:
        # write in whatever format necessary
        write_func = fmt_to_write_func[fmt]
        if fmt == 'pickle':
            output_dict = [v.__dict__ for v in serovar_predictions]
        else:
            output_dict = []
            for prediction in serovar_predictions:
                tmp = {}
                for k,v in prediction.__dict__.items():
                    if isinstance(v, (str, float, int)):
                        tmp[k] = v
                output_dict.append(tmp)
        write_func(fh, output_dict)
    except Exception as ex:
        logging.error(ex)
    finally:
        fh.close()

