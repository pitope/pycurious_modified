"""Continental scale area."""

import numpy as np
from pandas import DataFrame
from pycurious import CurieGrid, CurieOptimiseTanaka, CurieOptimiseBouligand, CurieOptimiseBansal
from pycurious.mapping import import_geotiff, export_geotiff
from areas.pietro_test.definition import runs, root_path, results_path

spacing = 50e3, 50e3
estimation_extent = (-2716000.0, 2785000.0, -2287000.0, 2423000.0) # x0, x1, y0, y1

if __name__ == '__main__':
    data_path = root_path / 'data/mag/continental/ADMAP2.0+mag4km_RCR_LCS1_infill_LCS1_LP_300km.tif'
    data, (x_0, x_1, y_0, y_1), projection, (d_x, d_y) = import_geotiff(data_path)
    curie = CurieGrid(data, x_0, x_1, y_0, y_1)
    xc_list, yc_list = curie.create_centroid_list(500e3, spacingX=spacing[0], spacingY=spacing[1],
                                                  subset=estimation_extent)

    for area_name, params in runs.items():
        data_path = root_path / 'data/mag/continental' / params['filename']
        data_path = data_path.resolve()
        if data_path.exists() is not True:
            raise FileNotFoundError

        data, (x_0, x_1, y_0, y_1), projection, (d_x, d_y) = import_geotiff(data_path)

        tanaka = CurieOptimiseTanaka(data, x_0, x_1, y_0, y_1)
        bouligand = CurieOptimiseBouligand(data, x_0, x_1, y_0, y_1)
        bansal = CurieOptimiseBansal(data, x_0, x_1, y_0, y_1)

        zt_range = params['zt_range']
        z0_range = params['z0_range']
        bansal_beta = params['bansal_beta']
        bouligand_k_range = params['bouligand_k_range']
        window_sizes = params['window_sizes']
        methods = params['methods']
        flight_height = params['flight_height']

        results = {
            'x': np.array(xc_list),
            'y': np.array(yc_list),
        }

        for window_size in window_sizes:
            w_fmt = f'{window_size / 1e3:g}' # in km
            corners = {
                f'{w_fmt}km_corner_0_x': xc_list - window_size / 2,
                f'{w_fmt}km_corner_0_y': yc_list - window_size / 2,
                f'{w_fmt}km_corner_1_x': xc_list - window_size / 2,
                f'{w_fmt}km_corner_1_y': yc_list + window_size / 2,
                f'{w_fmt}km_corner_2_x': xc_list + window_size / 2,
                f'{w_fmt}km_corner_2_y': yc_list + window_size / 2,
                f'{w_fmt}km_corner_3_x': xc_list + window_size / 2,
                f'{w_fmt}km_corner_3_y': yc_list - window_size / 2,
            }

            for method in methods:
                bouligand_run = False

                print(f'Running {method} on {len(xc_list)} centroids at {window_size} on {area_name}')
                products = {}
                if method == 'tanaka' or method == 'bansal':
                    if method == 'tanaka':
                        out = tanaka.optimise_routine(window_size, xc_list, yc_list,
                                                      zt_range, z0_range, taper=None, nan_fraction=0.5)
                    elif method == 'bansal':
                        out = bansal.optimise_routine(window_size, xc_list, yc_list,
                                                      zt_range, z0_range, beta=bansal_beta, taper=None, nan_fraction=0.5)

                    zt, z0, zt_int, z0_int, zt_stdev, z0_stdev = out
                    CPD, sigma_CPD = tanaka.calculate_CPD(zt, z0, zt_stdev, z0_stdev)

                    products = {
                        'z_top': -zt - flight_height,
                        'z_centroid': -z0 - flight_height,
                        'z_base': -CPD - flight_height,
                        'error_top': zt_stdev,
                        'error_centroid': z0_stdev,
                        'error_base': sigma_CPD,
                    }

                    for surface in 'top', 'base':
                        products[f'z_{surface}_min'] = products[f'z_{surface}'] - products[f'error_{surface}']
                        products[f'z_{surface}_max'] = products[f'z_{surface}'] + products[f'error_{surface}']

                elif method == 'bouligand':
                    if bouligand_run is True:
                        break
                    bouligand.reset_priors()
                    beta, zt, dz, C = bouligand.optimise_routine(window_size, xc_list, yc_list, taper=np.hanning, k_range=bouligand_k_range, nan_fraction=0.1)
                    products = {
                        'z_top': zt + flight_height,
                        'z_base': zt + dz + flight_height,
                        'beta': beta,
                        'C': C,
                    }
                    bouligand_run = True

                for key, value in products.items():
                    new_key = f'{w_fmt}km_{key}_{method}_{area_name}'
                    corners[new_key] = value

            results = results | corners

        df = DataFrame(results)
        df = df.round(2)
        df.to_csv(results_path / f'results_{area_name}.csv', index=False)