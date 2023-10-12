from pprint import pprint, pformat
from E_surface_intersection_finder import EsurfaceIntersectionFinder
import random
import logging
import sys
import vectors
import math
import re
import sympy as sp
from polynomial_formatter import format_polynomial
import argparse
import sys

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

def eval_oses(e_surf, ks):
    return [
        math.sqrt((sum(vectors.Vector(ks[i_loop_momentum])*sig for i_loop_momentum, sig in enumerate(
            os_e['loop_sig']) if sig != 0)+os_e['v_shift']).square()+os_e['m_squared'])
        for os_e in e_surf["onshell_propagators"]
    ]



def m_solve_intersection(e_surfaces, in_prec=128, out_prec=None,out='ms.in', ks=None):

    def format_for_m_solve(str_expr):
        # Below would be more precise but generate larger rationalized coefs 
        #sp.nsimplify(sp.parse_expr(sp_expr_str), rational=True, rational_conversion='exact').simplify().expand()
        # So better to expand first
        sp_expr=sp.nsimplify(sp.parse_expr(str_expr).expand(), rational=True, rational_conversion='exact').simplify()
        res = format_polynomial(str(sp_expr))
        return res


    n_loops = e_surfaces[0]['n_loops']
    n_E_surfaces = len(e_surfaces)

    coords=['x','y','z']
    ms_input_file = []
    ms_input_file.append(','.join( 
       ([','.join( 'k%d%s'%(i,j) for j in coords for i in range(1, n_loops+1)),] if ks is None else [])
        + [','.join( 'OSE%d%d'%(i,j) for i in range(1, n_E_surfaces+1) for j in range(1,len(e_surfaces[i-1]['onshell_propagators'])+1)),]
    ) )
    ms_input_file.append('0')
    ms_eqs = []
    for i_surf, e_surf in enumerate(e_surfaces):
        ms_eqs.append('+'.join('OSE%d%d'%(i_surf+1,j) for j in range(1,len(e_surf['onshell_propagators'])+1))+'%s'%format_for_m_solve(str(e_surf['E_shift'])))
        for i_prop, os_prop in enumerate(e_surf['onshell_propagators']):
            sp_expr_str = '%s+%s'%(
                '+'.join('(%s+(%s))**2'%( '+'.join(
                    '(%s%s%s%s)'%(
                        ('+' if s>0 else '-','k',i_sig+1,j) if ks is None else ('',ks[i_sig][i_coord],'','')
                    ) for i_sig, s in enumerate(os_prop['loop_sig']) if s!=0)
                ,
                os_prop['v_shift'][i_coord]
                ) for i_coord, j in enumerate(coords)),
                os_prop['m_squared']
            )
            ms_eqs.append(format_for_m_solve(sp_expr_str)+'-OSE%d%d^2'%(i_surf+1,i_prop+1))

    ms_input_file.append(',\n'.join(ms_eqs))

    with open(out,'w') as f:
        f.write('\n'.join(ms_input_file))

    def parse_ms_output(data):
        def replace_decimal(match):
            int_num = match.group(1)
            base = match.group(2)
            exp = match.group(3)
            return f"Decimal({int_num} / {base}**{exp})"
        modified_data = re.sub(r"(-?\d+) / (\d+)\^(\d+)", replace_decimal, data.replace(':',''))
        return eval(modified_data)

def random_three_two_loop_E_surfaces_examples(signatures, energy_scale=10., spatial_scale=1., mass_scale=1., seed=1, max_attempts=100, debug=False, use_m_solve=True):

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

        if use_m_solve:
            m_solve_intersection(E_surfaces,out='ms_candidate_case_%d.in'%n_attempts)
            logger.info("Wrote msolve input file to process case %d to file 'ms_candidate_case_%d.in'.\n> You can process it by running './msolve -p 128 -f ms_candidate_case_%d.in -o out.ms -v 2; python3 read_msolve.py'."%(n_attempts,n_attempts,n_attempts))

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
        logger.info("Found intersectin point = \n%s\nE-surf evaluations = %s" % (
            pformat([[
                ['%.16e' % k_i for k_i in k]
            ] for k in intersection_point]), ['%.16e' % eval_e_surf(e_surf, intersection_point) for e_surf in E_surfaces]))
        #res = m_solve_intersection(E_surfaces,ks=intersection_point, out='ms_test.in')
        return

    if n_attempts == 0:
        logger.warning(
            "Could not find a combination of random E-surfaces with an interior point after %d attempts" % max_attempts)
    else:
        logger.warning(
            "Could try to find interesection of E-surfaces %d times but never found one" % n_attempts)


if __name__ == '__main__':

    argparse = argparse.ArgumentParser("Test script for finding intersection of E-surfaces")
    argparse.add_argument("--test-case", help="Test case to run", default='1L_3E')
    args = argparse.parse_args()

    match args.test_case:
        case "2L_3E":
            signatures = [
                [[1, 0], [0, 1], [1, 1]],
                [[1, 0], [0, 1], [1, 1]],
                [[1, 0], [0, 1], [1, 1]]
            ]
        case "1L_3E":
            signatures = [
                [[1,], [1,]],
                [[1,], [1,]],
                [[1,], [1,]]
            ]
        case "2L_5E":
            signatures = [
                [[1, 0], [0, 1], [1, 1]],
                [[1, 0], [0, 1], [1, 1]],
                [[1, 0], [0, 1], [1, 1]],
                [[1, 0], [0, 1], [1, 1]],
                [[1, 0], [0, 1], [1, 1]],
            ]
        case _:
            logger.critical("Unknown test case %s" % args.test_case)
            sys.exit(1)
    
    DEBUG = logger.level == logging.DEBUG
    random_three_two_loop_E_surfaces_examples(
        signatures,
        seed=1,
        energy_scale=3., spatial_scale=1., mass_scale=1.,
        max_attempts=100, debug=DEBUG
    )
