from pprint import pprint, pformat
from E_surface_intersection_finder import EsurfaceIntersectionFinder
import random
import logging
import sys
import vectors
import math
logging.basicConfig(
    format='%(asctime)s :: %(levelname)s:%(message)s', level=logging.DEBUG)

logger = logging.getLogger('Test')

try:
    import cvxpy
    import scipy.optimize as optimize
except:
    logger.critical(
        "Error: could not import package cvxpy and scipy.optimize necessary for processing the E-surfaces for the IR profile. Make sure it is installed.")
    sys.exit(1)


def eval_e_surf(e_surf, ks):
    return sum(
        math.sqrt((sum(vectors.Vector(ks[i_loop_momentum])*sig for i_loop_momentum, sig in enumerate(
            os_e['loop_sig']) if sig != 0)+os_e['v_shift']).square()+os_e['m_squared'])
        for os_e in e_surf["onshell_propagators"]
    )+e_surf['E_shift']


def random_three_two_loop_E_surfaces_examples(signatures, energy_scale=10., spatial_scale=1., mass_scale=1., seed=1, max_attempts=100, debug=False):

    n_loops = len(signatures[0][0])
    n_E_surfaces = len(signatures)

    random.seed(seed)

    n_attempts = 0
    # Look for E-surfaces that likely have an intersection
    for i_attempt in range(max_attempts):
        E_surfaces = []
        cvxpy_coordinates = [cvxpy.Variable(3) for _ in range(n_loops)]

        for i_surf in range(n_E_surfaces):
            qs = [
                vectors.LorentzVector([
                    energy_scale*random.random(),
                    spatial_scale*2*(random.random()-0.5),
                    spatial_scale*2*(random.random()-0.5),
                    spatial_scale*2*(random.random()-0.5)
                ]) for i_p in range(len(signatures[i_surf]))
            ]
            masses = [mass_scale*random.random()
                      for i_p in range(len(signatures[i_surf]))]

            e_surf = {
                'id': i_surf,
                'E_surface_key': i_surf,
                'n_loops': n_loops,
                'onshell_propagators': [
                    {
                        'name': 'OS[%d][%d]' % (i_surf, i_p),
                        'square_root_sign': 1,
                        'energy_shift_sign': 1,
                        'loop_sig': tuple(sig),
                        'v_shift': q.space(),
                        'm_squared': m**2
                    } for i_p, (q, m, sig) in enumerate(zip(qs, masses, signatures[i_surf]))
                ],
                'E_shift': - sum(q[0] for q in qs)
            }

            e_surf['cxpy_expression'] = sum(
                cvxpy.norm(cvxpy.hstack([
                    math.sqrt(os_e['m_squared']),
                    sum(cvxpy_coordinates[i_loop_momentum] * sig for i_loop_momentum,
                        sig in enumerate(os_e['loop_sig']) if sig != 0)+list(os_e['v_shift'])
                ]), 2)
                for os_e in e_surf["onshell_propagators"]
            )+e_surf['E_shift']

            E_surfaces.append(e_surf)

        # First make sure there is an interior point to all E-surfaces
        p = cvxpy.Problem(
            cvxpy.Minimize(1),
            [E_surf['cxpy_expression'] <= 0 for E_surf in E_surfaces])

        try:
            cvxpyres = p.solve()
        except cvxpy.SolverError:
            p.solve(solver=cvxpy.SCS, eps=1e-9)
        except Exception as e:
            pass

        if p.status == cvxpy.OPTIMAL:
            # k_sol = [[float(c.value[0]), float(c.value[1]), float(
            #     c.value[2])] for c in cvxpy_coordinates]
            # logger.info("Interior point = %s\nE-surf evaluations = %s" % (
            #       k_sol, ['%.16e'%eval_e_surf(e_surf, k_sol) for e_surf in E_surfaces]))
            pass
        else:
            # Abort as if there is no interior point there, cannot be an intersection
            continue

        # Now find a point well in the interior by minimizing the sum of the E-surfaces evals
        p = cvxpy.Problem(
            cvxpy.Minimize(sum(E_surf['cxpy_expression']
                           for E_surf in E_surfaces)),
            [E_surf['cxpy_expression'] <= 0 for E_surf in E_surfaces])
        interior_point = [[float(c.value[0]), float(
            c.value[1]), float(c.value[2])] for c in cvxpy_coordinates]

        n_attempts += 1
        logger.info(
            "Attempting to find intersection for candidate case #%d" % n_attempts)
        logger.info("E-surface evaluation at interior point = \n%s" %
                    ['%.16e' % eval_e_surf(e_surf, interior_point) for e_surf in E_surfaces])

        intersection_point = EsurfaceIntersectionFinder(
            E_surfaces, cvxpy_coordinates, energy_scale, debug=debug).find_intersection()
        if intersection_point is None:
            logger.info(
                "Could not find an intersection for this candidate case #%d" % n_attempts)
            continue
        logger.info("Found intersectin point = %s\nE-surf evaluations = %s" % (
            intersection_point, ['%.16e' % eval_e_surf(e_surf, intersection_point) for e_surf in E_surfaces]))
        return

    if n_attempts == 0:
        logger.warning(
            "Could not find a combination of random E-surfaces with an interior point after %d attempts" % max_attempts)
    else:
        logger.warning(
            "Could try to find interesection of E-surfaces %d times but never found one" % n_attempts)


if __name__ == '__main__':

    DEBUG = logger.level == logging.DEBUG
    random_three_two_loop_E_surfaces_examples(
        [
            [[1, 0], [0, 1], [1, 1]],
            [[1, 0], [0, 1], [1, 1]],
            [[1, 0], [0, 1], [1, 1]]
        ],
        seed=1,
        energy_scale=3., spatial_scale=1., mass_scale=1.,
        max_attempts=100, debug=DEBUG
    )
