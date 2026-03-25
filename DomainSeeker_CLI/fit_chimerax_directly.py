#!/usr/lib/ucsf-chimerax-daily/bin/python3.9

import os, argparse
import numpy as np

import chimerax
from chimerax.core.session import Session, register_misc_commands
from chimerax.core import core_settings, attributes, startup, toolshed, version, tools, undo
from chimerax.core.commands import run

def parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('map_path')
    parser.add_argument('output_subdir')
    parser.add_argument('ref_map_threshold')
    parser.add_argument('resolution')
    parser.add_argument('n_search')
    parser.add_argument('domain_path', nargs="+")
    return parser.parse_args()


# domain_path=sys.argv[1]
# map_path=sys.argv[2]
# output_subdir=sys.argv[3]
# ref_map_threshold=float(sys.argv[4])
# resolution=float(sys.argv[5])
# n_search=int(sys.argv[6])


def session_setup() -> Session:
    from chimerax.core import version
    if version == '1.6.dev202302040126':
        import appdirs
        chimerax.app_dirs = ad = appdirs.AppDirs("ChimeraX", appauthor="",
                                                version="")
        chimerax.app_dirs_unversioned = adu = appdirs.AppDirs("ChimeraX", appauthor="")
        rootdir = "/"
        chimerax.app_bin_dir = os.path.join(rootdir, "bin")
        chimerax.app_data_dir = os.path.join(rootdir, "share")
        chimerax.app_lib_dir = os.path.join(rootdir, "lib")
        sess = Session("ChimeraX", silent=True)
        core_settings.init(sess)
        register_misc_commands(sess)
        attributes.RegAttrManager(sess)
        from chimerax.core.nogui import NoGuiLog
        sess.logger.add_log(NoGuiLog())
        startup.run_user_startup_scripts(sess)
        toolshed.init(
                sess.logger,
                debug=sess.debug,
                check_available=False,
                remote_url=None,
                session=sess
            )
        sess.toolshed = toolshed.get_toolshed()
        sess.toolshed.bootstrap_bundles(sess, False)
        sess.tools = tools.Tools(sess, first=True)
        sess.undo = undo.Undo(sess, first=True)
        return sess
    else:
        from chimerax.core.__main__ import _set_app_dirs, dedup_sys_path
        dedup_sys_path()
        _set_app_dirs(version)
        sess = Session(silent=True)
        core_settings.init(sess)
        register_misc_commands(sess)
        attributes.RegAttrManager(sess)
        from chimerax.core.nogui import NoGuiLog
        sess.logger.add_log(NoGuiLog())
        startup.run_user_startup_scripts(sess)
        toolshed.init(
                sess.logger,
                debug=sess.debug,
                check_available=False,
                remote_url=None,
                session=sess
            )
        sess.toolshed = toolshed.get_toolshed()
        sess.toolshed.bootstrap_bundles(sess, False)
        sess.tools = tools.Tools(sess, first=True)
        sess.undo = undo.Undo(sess, first=True)

        return sess



def fit_single_structure(session: Session, domain_path: str, output_subdir: str, ref_map_threshold: float, resolution: float, n_search: int):
    '''
    map has to be preloaded into session before passing session to fit_single_structure
    '''
    run(session,f'open {domain_path}')
    fits=run(session,f"fitmap #3 inmap #2 resolution {resolution} search {n_search}")
    # run(sess, 'open /home/ek/code/DomainSeeker/Example/em_data/A.mrc')
    # run(session, 'open /cds/ban/ekummerant/domainseeker/DomainSeeker/Example/em_data/A.mrc')
    # run(sess, 'open /home/ek/code/DomainSeeker/Example/pdb/Q6X6Z7.pdb')
    # run(session, "open /cds/ban/ekummerant/domainseeker/DomainSeeker/DomainSeeker_CLI/test/domains/A0A1Y7VMI0_D0.pdb")
    # fits = run(sess, 'fitmap #2 inmap #1 resolution 6.5 search 200')
    log_data=[]
    for fit in fits:
        transform_string_list=[]
        for element in fit.model_transforms()[0][1].matrix.reshape(-1):
            transform_string_list.append(f"{element:10.5f}")
        if fit.correlation() >= 0.00001:
            log_data.append([f"{fit.correlation():8.5f}",f"{fit.hits():8d}"]+transform_string_list)

    # write transformation info into log
    print(len(log_data))
    if len(log_data) >=1:
        os.makedirs(output_subdir,exist_ok=True)
        fitlog_subdir=os.path.join(output_subdir,"fitlogs")
        os.makedirs(fitlog_subdir,exist_ok=True)
        log_path=os.path.join(fitlog_subdir,os.path.basename(domain_path).replace('pdb','log'))
        print(f"{output_subdir=}")
        print(f"{log_path=}")
        np.savetxt(log_path,log_data,fmt="%s")
    
    run(session,f'close ~#1,2')



def main():
    argument = parser()

    sess = session_setup()
    run(sess,f'open {argument.map_path}')
    run(sess,f"volume threshold #1 minimum {float(argument.ref_map_threshold)} set 0")
    run(sess,f'volume #2 level {float(argument.ref_map_threshold)}')

    for domain in argument.domain_path:
        fit_single_structure(
            sess,
            domain,
            argument.output_subdir,
            float(argument.ref_map_threshold),
            float(argument.resolution),
            int(argument.n_search)
        )

if __name__ == "__main__":
    main()