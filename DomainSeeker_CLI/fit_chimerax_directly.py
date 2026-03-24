from chimerax.core.session import Session
from chimerax.core import core_settings
from chimerax.core.session import register_misc_commands
from chimerax.core import attributes
from chimerax.core.commands import run
from chimerax.core import startup
from chimerax.core import toolshed
from chimerax.core import version
from chimerax.core.__main__ import _set_app_dirs, dedup_sys_path
from chimerax.core import tools
from chimerax.core import undo

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
run(sess, 'open /home/ek/code/DomainSeeker/Example/em_data/A.mrc')
run(sess, 'open /home/ek/code/DomainSeeker/Example/pdb/Q6X6Z7.pdb')
run(sess, 'fitmap #2 inmap #1 resolution 6.5 search 200')
