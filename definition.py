from pathlib import Path
from cpd_utils import __path__ as paths
from copy import copy

root_path = Path(paths[0]) / '../..'

runs = {
    'Dziadek': {
        'filename': 'ADMAP2.0+mag4km_RCR_LCS1_infill_LCS1_LP_300km.tif',
        'flight_height': 4.0,
        'methods': ('tanaka',), #, 'bouligand'),
        'zt_range': (0.005, 0.014),
        'z0_range': (0.001, 0.005),
        'bouligand_k_range': (0.01, 0.5),
        'bansal_beta': 1.2,
        'window_sizes': (200e3, 300e3),
    }
}

"""
    'LCS1': {
        'filename': 'LCS-1_TMI_LP300km.tif',
        'flight_height': 0.0,
        'methods': ('tanaka','bansal'), #, 'bouligand'),
        'zt_range': (0.0044, 0.026),
        'z0_range': (0.0044, 0.026),
        'bouligand_k_range': (0.0, 0.1),
        'bansal_beta': 1.2,
        'window_sizes': (500e3, 1000e3),
    },
    'Martos-E': {
        'filename': 'ADMAP2.0+mag4km_RCR_LCS1_infill_LCS1_LP_300km.tif',
        'flight_height': 4.0,
        'methods': ('bansal', 'tanaka'), #, 'bouligand'),
        'zt_range': (0.041, 0.064),
        'z0_range': (0.0044, 0.026),
        'bouligand_k_range': (0.01, 0.5),
        'bansal_beta': 1.2,
        'window_sizes': (200e3, 300e3, 350e3),
    },
    'Martos-W': {
        'filename': 'ADMAP2.0+mag4km_RCR_LCS1_infill_LCS1_LP_300km.tif',
        'flight_height': 4.0,
        'methods': ('bansal', 'tanaka'), #, 'bouligand'),
        'zt_range': (0.051, 0.076),
        'z0_range': (0.01, 0.044),
        'bouligand_k_range': (0.01, 0.5),
        'bansal_beta': 1.2,
        'window_sizes': (200e3, 300e3, 350e3),
    },
    'Dziadek': {
        'filename': 'ADMAP2.0+mag4km_RCR_LCS1_infill_LCS1_LP_300km.tif',
        'flight_height': 4.0,
        'methods': ('bansal', 'tanaka'), #, 'bouligand'),
        'zt_range': (0.031, 0.088),
        'z0_range': (0.0063, 0.031),
        'bouligand_k_range': (0.01, 0.5),
        'bansal_beta': 1.2,
        'window_sizes': (200e3, 300e3),
    }
    
       'West_b3.0': {
        'filename': 'ADMAP2.0+mag4km_RCR_LCS1_infill_LCS1_LP_300km.tif',
        'flight_height': 4.0,
        'methods': ('bansal',),  # , 'bouligand'),
        'zt_range': (0.051, 0.076),
        'z0_range': (0.01, 0.044),
        'bouligand_k_range': (0.01, 0.5),
        'bansal_beta': 3.0,
        'window_sizes': (200e3,),
    }
}

runs['Martos-W'] = copy(runs['Martos-E'])
runs['Martos-W']['zt_range'] = [0.051, 0.076]
runs['Martos-W']['z0_range'] = [0.010, 0.045]"""


results_path = root_path / 'results/areas/pietro_test'