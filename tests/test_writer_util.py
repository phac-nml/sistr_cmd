from sistr_cmd.writers import to_dict
import numpy as np
import pandas as pd

class Test1:
    a = 1
    b = 3.2
    c = 'abc'
    z = False

    def __init__(self):
        self.i = 9
        self.j = 0.56
        self.k = 'test'
        self.l = [1, 2, 3]
        self.d = {'a':1, 'b':2}

class Test2:
    q = 6
    w = True
    e = 5.678
    r = [8.9, 10.1, 2.3, 4.5]

    def __init__(self, test1):
        self.test1 = test1
        a = 6
        b = 7.8
        c = 'string'
        d = {'key1': 'val1', 'key2': 'val2'}
        l = ['a', 'b', 'c']


def test_complex_obj_to_dict():
    t1 = Test1()
    t2 = Test2(t1)
    t2_dict = to_dict(t2, 0)
    assert isinstance(t2_dict, dict)
    exp_t2_dict = {'q': 6,
                   'test1': {'a': 1,
                             'c': 'abc',
                             'b': 3.2,
                             'd': {'a': 1,
                                   'b': 2},
                             'i': 9,
                             'k': 'test',
                             'j': 0.56,
                             'l': [1, 2, 3],
                             'z': False},
                   'r': [8.9, 10.1, 2.3, 4.5],
                   'e': 5.678,
                   'w': True}
    assert t2_dict == exp_t2_dict

    exp_t2_dict_depth2 = {'q': 6,
                   'test1': {'a': 1,
                             'c': 'abc',
                             'b': 3.2,
                             'i': 9,
                             'k': 'test',
                             'j': 0.56,
                             'z': False},
                   'r': [8.9, 10.1, 2.3, 4.5],
                   'e': 5.678,
                   'w': True}
    t2_dict_depth2 = to_dict(t2, 0, depth_threshold=2)
    assert t2_dict_depth2 == exp_t2_dict_depth2

    exp_t2_dict = {
                   'r': [8.9, 10.1, 2.3, 4.5],
                   'e': 5.678,
                   'w': True}
    t2_dict = to_dict(t2, 0, exclude_keys={'test1', 'q'})
    assert t2_dict == exp_t2_dict

    exp_t2_dict = {'q': 6,
                   'test1': {
                             'c': 'abc',
                             'b': 3.2,
                             'd': {
                                   'b': 2},
                             'i': 9,
                             'k': 'test',
                             'j': 0.56,
                             'l': [1, 2, 3],
                             'z': False},
                   'r': [8.9, 10.1, 2.3, 4.5],
                   'e': 5.678,
                   'w': True}
    t2_dict = to_dict(t2, 0, exclude_keys={'a'})
    assert t2_dict == exp_t2_dict


def test_numpy_numbers_to_dict():
    import json
    class ABC:
        def __init__(self, a, b, c):
            self.a = a
            self.b = b
            self.c = c

    t = [{'a': 1.90,
          'b': 7.8,
          'c': 99},
         {'a': 1.90,
          'b': 80.6,
          'c': 9999},
         {'a': 13.90,
          'b': 78.5,
          'c': 999},]
    df = pd.DataFrame(t)
    for i, r in df.iterrows():
        c = r.c
        c = c.astype(np.int64)
        abc = ABC(r.a, r.b, c)
        assert json.dumps(to_dict(abc, 0), sort_keys=True) == json.dumps(t[i], sort_keys=True)
