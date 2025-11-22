# 全局异常处理
import sys

def custom_excepthook(exc_type, exc_value, exc_traceback):
    """自定义异常钩子"""
    # 保存错误和额外信息到错误日志
    import os
    import traceback
    import datetime
    import platform
    log_str=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'\n'

    # 分割线
    log_str+='-'*80+'\n'

    # 系统信息
    log_str+=f"System: {platform.platform()}\n"
    # 解释器信息
    log_str+=f"Python version: {sys.version}\n"
    log_str+=f"Python executable: {sys.executable}\n"

    # 分割线
    log_str+='-'*80+'\n'

    # 工作目录
    log_str+=f"Working directory: {os.getcwd()}\n"
    # 调用命令
    log_str+=f"Command: {' '.join(sys.argv)}\n"

    # 分割线
    log_str+='-'*80+'\n'
    
    # 异常信息
    log_str+=f"Exception type: {exc_type.__name__}\n"
    log_str+=f"Exception value: {exc_value}\n"
    log_str+=f"Traceback:\n{''.join(traceback.format_tb(exc_traceback))}\n"
    # 保存日志
    error_log_path = sys.argv[1]
    with open(error_log_path, 'a') as f:
        f.write(log_str+'\n')
        # 结尾分割线
        f.write('='*80+'\n')
    # 退出程序
    sys.exit(1)

sys.excepthook = custom_excepthook
