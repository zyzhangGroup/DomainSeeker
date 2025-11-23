import os, sys
from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QTabWidget, QSpacerItem,
                          QLabel, 
                          QLineEdit, QTextEdit, 
                          QPushButton, 
                          QFileDialog, 
                          QFrame, QSizePolicy,
                          QScrollArea )
from Qt.QtCore import Qt
import numpy as np
from chimerax.core.commands import run
from chimerax.geometry.place import Place
import datetime
import subprocess

script_dir=os.path.dirname(os.path.realpath(__file__))

# 不同系统的ChimeraX
if sys.platform == 'win32':
    ChimeraX="ChimeraX-console"     # windows
else:
    ChimeraX="chimerax"     # Linux

class DomainSeeker(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = False         # We do not save/restore in sessions temporarily
    # No help document temporarily
    # help = "help:user/tools/DomainSeeker.html"

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string
        
        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "DomainSeeker"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        #
        # Note that by default, tool windows are only hidden rather than
        # destroyed when the user clicks the window's close button.  To change
        # this behavior, specify 'close_destroys=True' in the MainToolWindow
        # constructor.
        self.tool_window = MainToolWindow(self)

        # global options
        self.prior_results_loaded = False
        self.posterior_results_loaded = False

        # Our user interface is simple enough that we could probably inline
        # the code right here, but for any kind of even moderately complex
        # interface, it is probably better to put the code in a method so
        # that this __init__ method remains readable.
        self._build_ui()


    def _build_ui(self):
        # Put our widgets in the tool window

        # The main layout contains other layouts for seperate parts of 
        # our tool interface.
        main_layout = QVBoxLayout()

        # 添加标签页控件
        tab_widget = QTabWidget()
        main_layout.addWidget(tab_widget)

        # 添加计算页
        computation_tab = self._create_computation_tab()
        tab_widget.addTab(computation_tab, "Computation")

        # 添加展示页
        presentation_tab = self._create_presentation_tab()
        tab_widget.addTab(presentation_tab, "Presentation")
        
        # Set the main layout as the contents of our window
        self.tool_window.ui_area.setLayout(main_layout)

        # Show the window on the user-preferred side of the ChimeraX
        # main window
        self.tool_window.manage('side')

    # 创建计算标签页
    def _create_computation_tab(self):
        # Create a new QWidget to hold our widgets
        calculation_tab = QWidget()
        calculation_tab_base_layout = QVBoxLayout(calculation_tab)

        # 创建滚动区域
        computation_scroll_area = QScrollArea()
        computation_scroll_area.setObjectName("computation_scroll_area")
        computation_scroll_area.setWidgetResizable(True)
        # 添加到计算标签页布局
        calculation_tab_base_layout.addWidget(computation_scroll_area)

        # 创建滚动区域容器
        computation_scroll_area_container = QWidget()
        computation_scroll_area_container.setObjectName("result_scroll_area_container")
        # 垂直方向尽可能大
        computation_scroll_area_container.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        computation_scroll_area.setWidget(computation_scroll_area_container)

        # 创建网格布局，并添加到滚动区域容器
        calculation_layout = QVBoxLayout(computation_scroll_area_container)

        global_options_layout = QVBoxLayout()
        file_fetching_layout = QVBoxLayout()
        domain_parsing_layout = QVBoxLayout()
        fitting_scoring_layout = QVBoxLayout()
        prior_probability_layout = QVBoxLayout()
        posterior_probability_layout = QVBoxLayout()

        # --------------------------------------------------------------------------------------------------------

        # Global options layout

        global_options_head = QLabel("Global options")

        global_options_layout.addWidget(global_options_head)

        global_options_body = QGridLayout()

        # Project directory
        global_options_body.addWidget(QLabel("Project directory:"), 0, 0)

        self.project_directory_text = QLineEdit()
        # 禁止输入，只能通过按钮选择，因为需要初始化一些信息
        self.project_directory_text.setReadOnly(True)
        global_options_body.addWidget(self.project_directory_text, 0, 1)

        project_directory_select_button = QPushButton("Select Directory")
        global_options_body.addWidget(project_directory_select_button, 0, 2)
        project_directory_select_button.clicked.connect(lambda: self._initialize_project(self.project_directory_text))

        # Map directory
        global_options_body.addWidget(QLabel("Map directory:"), 1, 0)

        self.map_directory_text = QLineEdit()
        global_options_body.addWidget(self.map_directory_text, 1, 1)

        map_directory_select_button = QPushButton("Select Directory")
        global_options_body.addWidget(map_directory_select_button, 1, 2)
        map_directory_select_button.clicked.connect(lambda: self._select_directory(self.map_directory_text,
                                                                                   self.project_directory_text.text()))
        
        # Pdb directory
        global_options_body.addWidget(QLabel("Pdb directory:"), 2, 0)

        self.pdb_directory_text = QLineEdit()
        global_options_body.addWidget(self.pdb_directory_text, 2, 1)

        pdb_directory_select_button = QPushButton("Select Directory")
        global_options_body.addWidget(pdb_directory_select_button, 2, 2)
        pdb_directory_select_button.clicked.connect(lambda: self._select_directory(self.pdb_directory_text,
                                                                                   self.project_directory_text.text()))
        
        # Pae directory
        global_options_body.addWidget(QLabel("Pae directory:"), 3, 0)

        self.pae_directory_text = QLineEdit()        
        global_options_body.addWidget(self.pae_directory_text, 3, 1)

        pae_directory_select_button = QPushButton("Select Directory")
        global_options_body.addWidget(pae_directory_select_button, 3, 2)
        pae_directory_select_button.clicked.connect(lambda: self._select_directory(self.pae_directory_text,
                                                                                   self.project_directory_text.text()))
        
        
        # Domain directory
        global_options_body.addWidget(QLabel("Domain directory:"), 4, 0)

        self.domain_directory_text = QLineEdit()
        global_options_body.addWidget(self.domain_directory_text, 4, 1)

        domain_directory_select_button = QPushButton("Select Directory")
        global_options_body.addWidget(domain_directory_select_button, 4, 2)
        domain_directory_select_button.clicked.connect(lambda: self._select_directory(self.domain_directory_text,
                                                                                      self.project_directory_text.text()))
        
        # Fitout directory
        global_options_body.addWidget(QLabel("Fitout directory:"), 5, 0)

        self.fitout_directory_text = QLineEdit()
        global_options_body.addWidget(self.fitout_directory_text, 5, 1)

        fitout_directory_select_button = QPushButton("Select Directory")
        global_options_body.addWidget(fitout_directory_select_button, 5, 2)
        fitout_directory_select_button.clicked.connect(lambda: self._select_directory(self.fitout_directory_text,
                                                                                      self.project_directory_text.text()))
        
        global_options_layout.addLayout(global_options_body)

        # --------------------------------------------------------------------------------------------------------

        # File fetching layout

        file_fetching_head = QLabel("Fetch pdb and pae files from AFDB")

        file_fetching_layout.addWidget(file_fetching_head)


        file_fetching_body = QHBoxLayout()

        file_fetching_body.addWidget(QLabel("Fetch proteins in file:"))

        protein_list_file_path_text = QLineEdit()
        file_fetching_body.addWidget(protein_list_file_path_text)

        protein_list_file_select_button = QPushButton("Select File")
        file_fetching_body.addWidget(protein_list_file_select_button)
        protein_list_file_select_button.clicked.connect(lambda: self._select_file(protein_list_file_path_text, 
                                                                                  self.project_directory_text.text()))

        fetch_file_button = QPushButton("Fetch files")
        file_fetching_body.addWidget(fetch_file_button)
        fetch_file_button.clicked.connect(lambda: self._fetch_pdb_and_pae_files(protein_list_file_path_text.text(), 
                                                                                self.project_directory_text.text(),
                                                                                self.pdb_directory_text.text(),
                                                                                self.pae_directory_text.text()))
        file_fetching_layout.addLayout(file_fetching_body)

        # --------------------------------------------------------------------------------------------------------

        # Domain parsing layout

        domain_parsing_head = QLabel("Parse proteins into domains based on PAE")

        domain_parsing_layout.addWidget(domain_parsing_head)


        domain_parsing_body = QGridLayout()

        # plddt_cutoff
        domain_parsing_body.addWidget(QLabel("plddt_cutoff"), 0, 0)
        plddt_cutoff_text = QLineEdit()
        plddt_cutoff_text.setText("70")
        domain_parsing_body.addWidget(plddt_cutoff_text, 1, 0)

        # pae_cutoff
        domain_parsing_body.addWidget(QLabel("pae_cutoff"), 0, 1)
        pae_cutoff_text = QLineEdit()
        pae_cutoff_text.setText("5")
        domain_parsing_body.addWidget(pae_cutoff_text, 1, 1)

        # clique_cutoff
        domain_parsing_body.addWidget(QLabel("clique_cutoff"), 0, 2)
        clique_cutoff_text = QLineEdit()
        clique_cutoff_text.setText("4")
        domain_parsing_body.addWidget(clique_cutoff_text, 1, 2)

        # min_dege_ratio
        domain_parsing_body.addWidget(QLabel("min_dege_ratio"), 2, 0)
        min_dege_ratio_text = QLineEdit()
        min_dege_ratio_text.setText("0.6")
        domain_parsing_body.addWidget(min_dege_ratio_text, 3, 0)

        # min_common_nodes_ratio
        domain_parsing_body.addWidget(QLabel("min_common_nodes_ratio"), 2, 1)
        min_common_nodes_ratio_text = QLineEdit()
        min_common_nodes_ratio_text.setText("0.5")
        domain_parsing_body.addWidget(min_common_nodes_ratio_text, 3, 1)

        # min_domain_size
        domain_parsing_body.addWidget(QLabel("min_domain_size"), 2, 2)
        min_domain_size_text = QLineEdit()
        min_domain_size_text.setText("40")
        domain_parsing_body.addWidget(min_domain_size_text, 3, 2)

        # n_process
        domain_parsing_body.addWidget(QLabel("n_process"), 0, 3)
        domain_parsing_n_process_text = QLineEdit()
        domain_parsing_n_process_text.setText("1")
        domain_parsing_body.addWidget(domain_parsing_n_process_text, 1, 3)

        # Submit button
        domain_parsing_button = QPushButton("Parse domains")
        domain_parsing_body.addWidget(domain_parsing_button,3,3)
        domain_parsing_button.clicked.connect(lambda: self._parse_domains(self.project_directory_text.text(),
                                                                          self.pdb_directory_text.text(),
                                                                          self.pae_directory_text.text(),
                                                                          self.domain_directory_text.text(),
                                                                          plddt_cutoff_text.text(),
                                                                          pae_cutoff_text.text(),
                                                                          clique_cutoff_text.text(),
                                                                          min_dege_ratio_text.text(),
                                                                          min_common_nodes_ratio_text.text(),
                                                                          min_domain_size_text.text(),
                                                                          domain_parsing_n_process_text.text()))

        domain_parsing_layout.addLayout(domain_parsing_body)

        # --------------------------------------------------------------------------------------------------------

        # Fitting and scoring layout

        fitting_scoring_head = QLabel("Fit domains into densities and score fitted domains")
        fitting_scoring_layout.addWidget(fitting_scoring_head)

        fitting_scoring_body = QGridLayout()

        # threshold
        fitting_scoring_body.addWidget(QLabel("threshold"), 0, 0)
        threshold_text = QLineEdit()
        threshold_text.setText("0.0")
        fitting_scoring_body.addWidget(threshold_text, 1, 0)

        # resolution
        fitting_scoring_body.addWidget(QLabel("resolution"), 0, 1)
        resolution_text = QLineEdit()
        resolution_text.setText("6.0")
        fitting_scoring_body.addWidget(resolution_text, 1, 1)

        # n_search
        fitting_scoring_body.addWidget(QLabel("n_search"), 0, 2)
        n_search_text = QLineEdit()
        n_search_text.setText("200")
        fitting_scoring_body.addWidget(n_search_text, 1, 2)

        # n_process
        fitting_scoring_body.addWidget(QLabel("n_process"), 0, 3)
        fitting_n_process_text = QLineEdit()
        fitting_n_process_text.setText("1")
        fitting_scoring_body.addWidget(fitting_n_process_text, 1, 3)

        # negtive_laplacian_cutoff
        fitting_scoring_body.addWidget(QLabel("negtive_laplacian_cutoff"), 2, 0)
        negtive_laplacian_cutoff_text = QLineEdit()
        negtive_laplacian_cutoff_text.setText("-0.001")
        fitting_scoring_body.addWidget(negtive_laplacian_cutoff_text, 3, 0)

        # positive_laplacian_cutoff
        fitting_scoring_body.addWidget(QLabel("positive_laplacian_cutoff"), 2, 1)
        positive_laplacian_cutoff_text = QLineEdit()
        positive_laplacian_cutoff_text.setText("0.001")
        fitting_scoring_body.addWidget(positive_laplacian_cutoff_text, 3, 1)


        # Submit button
        fitting_scoring_button = QPushButton("Fit & score")
        fitting_scoring_body.addWidget(fitting_scoring_button,3,3)
        fitting_scoring_button.clicked.connect(lambda: self._fit_and_score(self.project_directory_text.text(),
                                                                           self.map_directory_text.text(),
                                                                           threshold_text.text(),
                                                                           resolution_text.text(),
                                                                           n_search_text.text(),
                                                                           negtive_laplacian_cutoff_text.text(),
                                                                           positive_laplacian_cutoff_text.text(),
                                                                           fitting_n_process_text.text(),
                                                                           self.domain_directory_text.text(),
                                                                           self.fitout_directory_text.text()))
        
        fitting_scoring_layout.addLayout(fitting_scoring_body)

        # --------------------------------------------------------------------------------------------------------

        # Prior probability layout

        prior_probability_head = QLabel("Calculate prior probability of each fitted domain")
        prior_probability_layout.addWidget(prior_probability_head)

        prior_probability_body = QGridLayout()

        # box_num
        prior_probability_body.addWidget(QLabel("box_num"), 0, 0)
        box_num_text = QLineEdit()
        box_num_text.setText("10")
        prior_probability_body.addWidget(box_num_text, 1, 0)

        # min_data_per_box
        prior_probability_body.addWidget(QLabel("min_data_per_box"), 0, 1)
        min_data_per_box_text = QLineEdit()
        min_data_per_box_text.setText("50")
        prior_probability_body.addWidget(min_data_per_box_text, 1, 1)

        # relative_density_cutoff
        prior_probability_body.addWidget(QLabel("relative_density_cutoff"), 0, 2)
        relative_density_cutoff_text = QLineEdit()
        relative_density_cutoff_text.setText("0.01")
        prior_probability_body.addWidget(relative_density_cutoff_text, 1, 2)

        # zScore_offset
        prior_probability_body.addWidget(QLabel("zScore_offset"), 2, 0)
        zScore_offset_text = QLineEdit()
        zScore_offset_text.setText("15")
        prior_probability_body.addWidget(zScore_offset_text, 3, 0)


        # Submit button
        prior_probability_button = QPushButton("Calculate prior probability")
        prior_probability_body.addWidget(prior_probability_button,3,3)
        prior_probability_button.clicked.connect(lambda: self._calculate_prior_probability(self.project_directory_text.text(),
                                                                                           self.map_directory_text.text(),
                                                                                           self.fitout_directory_text.text(),
                                                                                           box_num_text.text(),
                                                                                           min_data_per_box_text.text(),
                                                                                           relative_density_cutoff_text.text(),
                                                                                           zScore_offset_text.text()))

        prior_probability_layout.addLayout(prior_probability_body)

        # --------------------------------------------------------------------------------------------------------

        # Posterior probability layout

        posterior_probability_head = QLabel("Integrate extra experimental data")

        posterior_probability_layout.addWidget(posterior_probability_head)

        posterior_probability_body = QVBoxLayout()

        # XL-MS data
        crosslink_layout = QVBoxLayout()

        crosslink_head = QLabel("XL-MS data")

        crosslink_layout.addWidget(crosslink_head)

        crosslink_body = QHBoxLayout()

        # Manege crosslink files, adding, presenting, and removing files
        crosslink_files_manager = QVBoxLayout()
        # 简单实现，用按钮选择多个文件，在多行文本框分行显示
        
        # A button to add crosslink files
        select_crosslink_files_button = QPushButton("Select XL-MS files")
        crosslink_files_manager.addWidget(select_crosslink_files_button)
        select_crosslink_files_button.clicked.connect(lambda: self._select_files(crosslink_files_text,
                                                                                 self.project_directory_text.text()))

        # A text box to show crosslink files
        crosslink_files_text = QTextEdit()
        crosslink_files_text.setReadOnly(True)
        crosslink_files_manager.addWidget(crosslink_files_text)

        # Add the crosslink files manager to the crosslink layout
        crosslink_body.addLayout(crosslink_files_manager,stretch=2)

        # Options for integrating XL-MS data
        crosslink_options_layout = QGridLayout()

        # threshold, acceptor_cutoff, donor_cutoff, evidence_strenght
        crosslink_options_layout.addWidget(QLabel("threshold"), 0, 0)
        post_threshold_text = QLineEdit()
        post_threshold_text.setText("0.0")
        crosslink_options_layout.addWidget(post_threshold_text, 0, 1)

        crosslink_options_layout.addWidget(QLabel("acceptor_cutoff"), 1, 0)
        acceptor_cutoff_text = QLineEdit()
        acceptor_cutoff_text.setText("0.00001")
        crosslink_options_layout.addWidget(acceptor_cutoff_text, 1, 1)

        crosslink_options_layout.addWidget(QLabel("donor_cutoff"), 2, 0)
        donor_cutoff_text = QLineEdit()
        donor_cutoff_text.setText("0.01")
        crosslink_options_layout.addWidget(donor_cutoff_text, 2, 1)

        crosslink_options_layout.addWidget(QLabel("evidence_strength"), 3, 0)
        evidence_strength_text = QLineEdit()
        evidence_strength_text.setText("10")
        crosslink_options_layout.addWidget(evidence_strength_text, 3, 1)


        # Add the crosslink options to the crosslink layout
        crosslink_body.addLayout(crosslink_options_layout,stretch=1)

        crosslink_layout.addLayout(crosslink_body)

        # Add the crosslink layout to the posterior probability body
        posterior_probability_body.addLayout(crosslink_layout)

        # Submit button
        posterior_probability_button = QPushButton("Integrate experimental data")
        posterior_probability_body.addWidget(posterior_probability_button)
        posterior_probability_button.clicked.connect(lambda: self._calculate_posterior_probability(self.project_directory_text.text(),
                                                                                                   self.domain_directory_text.text(),
                                                                                                   self.map_directory_text.text(),
                                                                                                   post_threshold_text.text(),
                                                                                                   self.fitout_directory_text.text(),
                                                                                                   acceptor_cutoff_text.text(),
                                                                                                   donor_cutoff_text.text(),
                                                                                                   evidence_strength_text.text(),
                                                                                                   crosslink_files_text.toPlainText().split("\n")))
        
        posterior_probability_layout.addLayout(posterior_probability_body)


        # --------------------------------------------------------------------------------------------------------

        # 底部弹簧spacer
        vertival_bottom_spacer_of_calculation_tab = QSpacerItem(0, 0, QSizePolicy.Minimum, QSizePolicy.Expanding)

        # --------------------------------------------------------------------------------------------------------

        # add all layouts to main layout

        calculation_layout.addLayout(global_options_layout)

        calculation_layout.addWidget(self.create_horizontal_line())

        calculation_layout.addLayout(file_fetching_layout)

        calculation_layout.addWidget(self.create_horizontal_line())

        calculation_layout.addLayout(domain_parsing_layout)

        calculation_layout.addWidget(self.create_horizontal_line())

        calculation_layout.addLayout(fitting_scoring_layout)

        calculation_layout.addWidget(self.create_horizontal_line())

        calculation_layout.addLayout(prior_probability_layout)

        calculation_layout.addWidget(self.create_horizontal_line())

        calculation_layout.addLayout(posterior_probability_layout)

        calculation_layout.addWidget(self.create_horizontal_line(thickness=4))

        calculation_layout.addItem(vertival_bottom_spacer_of_calculation_tab)

        return calculation_tab
    
    # project 初始化
    def _initialize_project(self, project_directory_text):
        # 选择目录
        self._select_directory(project_directory_text)
        # 检查project目录下是否存在error.log文件。如果不存在，创建一个，并记录日志创建时间
        error_log_path = os.path.join(project_directory_text.text(), "error.log")
        error_log_path = error_log_path.replace("\\", "/")
        if not os.path.exists(error_log_path):
            with open(error_log_path, "w") as f:
                f.write(f"Log created at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                # 开头分隔符
                f.write('='*80+'\n')
        # 将error文件路径保存到self中
        self.error_log_path = error_log_path
        
                
        

    # 添加展示页
    def _create_presentation_tab(self):
        # 创建展示页容器
        presentation_tab = QWidget()
        # 添加布局
        presentation_layout = QVBoxLayout(presentation_tab)
        # 顶部对齐
        presentation_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        result_presentation_layout = QVBoxLayout()

        # Result presentation layout

        # head层包含标题
        result_presentation_head = QHBoxLayout()
        result_presentation_head.addWidget(QLabel("Result presentation"))

        result_presentation_layout.addLayout(result_presentation_head)

        # body层包含主要部分
        result_presentation_body = QVBoxLayout()

        # 添加控制按钮层，grid布局
        result_presentation_control_layout = QGridLayout()

        # 初始化按钮
        result_initialization_button = QPushButton("Initialize results")
        result_presentation_control_layout.addWidget(result_initialization_button, 0, 0)
        result_initialization_button.clicked.connect(lambda: self._initialize_results(result_presentation_body,
                                                                                      self.map_directory_text.text(),
                                                                                      self.project_directory_text.text(),
                                                                                      self.domain_directory_text.text(),
                                                                                      self.fitout_directory_text.text()))
        
        # a button to get posterior results
        get_posterior_results_button = QPushButton("Get posterior results")
        result_presentation_control_layout.addWidget(get_posterior_results_button, 0, 1)
        get_posterior_results_button.clicked.connect(lambda: self._get_posterior_results(self.project_directory_text.text(),
                                                                                         self.fitout_directory_text.text()))

        # a button to update results by prior ranks
        update_results_by_prior_ranks_button = QPushButton("Update by prior ranks")
        result_presentation_control_layout.addWidget(update_results_by_prior_ranks_button, 0, 2)
        update_results_by_prior_ranks_button.clicked.connect(lambda: self._update_results_by_prior_ranks(self.project_directory_text.text(),
                                                                                                         self.domain_directory_text.text(),
                                                                                                         self.fitout_directory_text.text()))

        # a button to update results by posterior ranks
        update_results_by_posterior_ranks_button = QPushButton("Update by posterior ranks")
        result_presentation_control_layout.addWidget(update_results_by_posterior_ranks_button, 0, 3)
        update_results_by_posterior_ranks_button.clicked.connect(lambda: self._update_results_by_posterior_ranks(self.project_directory_text.text(),
                                                                                                                 self.domain_directory_text.text(),
                                                                                                                 self.fitout_directory_text.text()))
        
        # 将控制层添加到body
        result_presentation_body.addLayout(result_presentation_control_layout)

        # Add the result_presentation_body to the result presentation layout
        result_presentation_layout.addLayout(result_presentation_body)

        # Add the result_presentation_layout to the presentation layout
        presentation_layout.addLayout(result_presentation_layout)

        return presentation_tab


    def create_horizontal_line(self, style="solid", color="#cccccc", thickness=2):
        """Create a custom horizontal separator line"""
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        
        styles = {
            "solid": f"border: {thickness}px {style} {color};",
            "dashed": f"border: {thickness}px dashed {color};",
            "dotted": f"border: {thickness}px dotted {color};",
            "double": f"border: {thickness}px double {color};",
            }
        
        line.setStyleSheet(styles.get(style, styles["solid"]))
        return line

    def _select_directory(self, target_text_edit, start_directory=""):
        try:
            # If the start directory does not exist, reset it to an empty string
            if not os.path.exists(start_directory):
                start_directory = ""
            directory_path = QFileDialog.getExistingDirectory(None, 
                                                              "Select the project directory", 
                                                              start_directory)
            if directory_path:
                target_text_edit.setText(directory_path)
        except Exception as e:
            self.session.logger.error(f"Error during directory selection: {e}")

    def _select_file(self, target_text_edit, start_directory=""):
        try:
            # If the start directory does not exist, reset it to an empty string
            if not os.path.exists(start_directory):
                start_directory = ""
            file_path, _ = QFileDialog.getOpenFileName(None, 
                                                       "Select a file containing candidate proteins", 
                                                       start_directory, 
                                                       "Text files (*.txt *.csv)")
            if file_path:
                target_text_edit.setText(file_path)
        except Exception as e:
            self.session.logger.error(f"Error during file selection: {e}")

    def _select_files(self, target_text_edit, start_directory=""):
        try:
            # If the start directory does not exist, reset it to an empty string
            if not os.path.exists(start_directory):
                start_directory = ""
            file_paths, _ = QFileDialog.getOpenFileNames(None, 
                                                         "Select XL-MS files", 
                                                         start_directory, 
                                                         "Text files (*.txt *.csv)")
            if file_paths:
                target_text_edit.setText("\n".join(file_paths))
        except Exception as e:
            self.session.logger.error(f"Error during file selection: {e}")
    
    def _fetch_pdb_and_pae_files(self, protein_list_file_path, project_directory, pdb_directory = "", pae_directory = ""):
        # 检查protein_list_file_path
        if not protein_list_file_path:
            self.session.logger.error("Please select a file containing candidate proteins")
            return
        elif not os.path.exists(protein_list_file_path):
            self.session.logger.error(f"Protein list file {protein_list_file_path} does not exist")
            return
        # Check if the project directory is valid, if not, raise an error
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # Check if the pdb directory is valid
        if not pdb_directory:
            pdb_directory = os.path.join(project_directory, "pdb_files")
        elif not os.path.exists(pdb_directory):
            self.session.logger.error(f"PDB directory {pdb_directory} does not exist")
            return
        pdb_directory = pdb_directory.replace("\\", "/")
        # Check if the pae directory is valid
        if not pae_directory:
            pae_directory = os.path.join(project_directory, "pae_files")
        elif not os.path.exists(pae_directory):
            self.session.logger.error(f"PAE directory {pae_directory} does not exist")
            return
        pae_directory = pae_directory.replace("\\", "/")
        # 运行子进程
        self.session.logger.info("start fetching pdb and pae files")
        arg_list=[f'{script_dir}/fetch_pdb_pae.py',
                  self.error_log_path,
                  protein_list_file_path,
                  project_directory,
                  pdb_directory,
                  pae_directory]
        self.run_detatched_subprocess(arg_list)

    def _initialize_results(self, target_layout, map_directory, project_directory, domain_directory = "", fitout_dir = ""):
        # 生成空白结果框架
        result_layout = self._generate_blank_results(target_layout,map_directory)
        # 将result_layout添加到self中，以便后续更新
        self.result_layout = result_layout
        # 打开密度文件
        self._open_density_files(map_directory)
        # 初始化fitted_domain_models记录
        self.fitted_domain_models = [(None,None) for _ in range(len(self.density_names))] # (state_id, atomic_model)
        # 读取先验概率结果， 并初始化选择器
        self._get_prior_results(project_directory,fitout_dir)
        # 显示初始状态
        self._update_current_states()
        # 显示初始状态的先验结果
        self._update_prior_results()
        # 更新原子模型
        self._update_fitted_domains(project_directory,domain_directory,fitout_dir)
        

        
    # 建立空白结果grid
    def _generate_blank_results(self, target_layout,map_directory):
        # 创建滚动区域
        result_scroll_area = QScrollArea()
        result_scroll_area.setObjectName("result_scroll_area")
        result_scroll_area.setWidgetResizable(True)

        # 创建滚动区域容器
        result_scroll_area_container = QWidget()
        result_scroll_area_container.setObjectName("result_scroll_area_container")
        # 垂直方向尽可能大
        result_scroll_area_container.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        result_scroll_area.setWidget(result_scroll_area_container)

        # 创建网格布局，并添加到滚动区域容器
        result_layout = QGridLayout(result_scroll_area_container)
        result_layout.setObjectName("result_layout")
        # 顶端对齐
        result_layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # density, domain, fit, prior_prob,prior_rank, posterior_prob, posterior_rank
        result_layout.addWidget(QLabel("density"), 0, 0)
        result_layout.setColumnStretch(0, 2)
        result_layout.addWidget(QLabel("domain"), 0, 1)
        result_layout.setColumnStretch(1, 2)
        result_layout.addWidget(QLabel("fit"), 0, 2)
        result_layout.setColumnStretch(2, 2)
        result_layout.addWidget(QLabel("pri_prob"), 0, 3)
        result_layout.setColumnStretch(3, 2)
        result_layout.addWidget(QLabel("pri_rank"), 0, 4)
        result_layout.setColumnStretch(4, 2)
        result_layout.addWidget(QLabel("pos_prob"), 0, 5)
        result_layout.setColumnStretch(5, 2)
        result_layout.addWidget(QLabel("pos_rank"), 0, 6)
        result_layout.setColumnStretch(6, 2)
        
        # # 总体状态结果
        # # Add a label to show the aggregated result
        # result_layout.addWidget(QLabel("Whole"), 1, 0)

        # # prior_prob,prior_rank, posterior_prob, posterior_rank
        # aggregated_prior_prob_text = QLineEdit()
        # aggregated_prior_prob_text.setReadOnly(True)
        # result_layout.addWidget(aggregated_prior_prob_text, 1, 3)

        # aggregated_prior_rank_text = QLineEdit()
        # aggregated_prior_rank_text.setReadOnly(False)
        # result_layout.addWidget(aggregated_prior_rank_text, 1, 4)

        # aggregated_posterior_prob_text = QLineEdit()
        # aggregated_posterior_prob_text.setReadOnly(True)
        # result_layout.addWidget(aggregated_posterior_prob_text, 1, 5)

        # aggregated_posterior_rank_text = QLineEdit()
        # aggregated_posterior_rank_text.setReadOnly(False)
        # result_layout.addWidget(aggregated_posterior_rank_text, 1, 6)


        # 单密度结果
        # 获取所有density文件名
        self.density_names = [file_name[:-4] for file_name in os.listdir(map_directory) if file_name.endswith(".mrc")]

        for i, density_name in enumerate(self.density_names):
            result_layout.addWidget(QLabel(density_name), 2+i, 0)

            domain_name_text = QLineEdit()
            domain_name_text.setReadOnly(True)
            result_layout.addWidget(domain_name_text, 2+i, 1)

            fit_id_text = QLineEdit()
            fit_id_text.setReadOnly(True)
            result_layout.addWidget(fit_id_text, 2+i, 2)

            prior_prob_text = QLineEdit()
            prior_prob_text.setReadOnly(True)
            result_layout.addWidget(prior_prob_text, 2+i, 3)

            prior_rank_text = QLineEdit()
            prior_rank_text.setReadOnly(False)
            result_layout.addWidget(prior_rank_text, 2+i, 4)

            posterior_prob_text = QLineEdit()
            posterior_prob_text.setReadOnly(True)
            result_layout.addWidget(posterior_prob_text, 2+i, 5)

            posterior_rank_text = QLineEdit()
            posterior_rank_text.setReadOnly(False)
            result_layout.addWidget(posterior_rank_text, 2+i, 6)
        
        # 设置相对宽度
        result_layout.setColumnStretch(0, 2)
        result_layout.setColumnStretch(1, 2)
        result_layout.setColumnStretch(2, 2)
        result_layout.setColumnStretch(3, 2)
        result_layout.setColumnStretch(4, 2)
        result_layout.setColumnStretch(5, 2)
        result_layout.setColumnStretch(6, 2)
        
        # 将滚动区域添加到结果布局
        target_layout.addWidget(result_scroll_area)

        return result_layout


    # 获取先验概率结果
    def _get_prior_results(self, project_directory, fitout_dir = ""):
        # 判断fitout_dir是否为空，为空则设置为默认值
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        # 验证fitout_dir是否存在
        if not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        # 记录每个密度的状态
        states_of_densities=[]
        states_to_id_dict=[]
        # 先验概率列表
        prior_prob_list = []
        # 记录每个状态的排名
        prior_sorted_state_ids = []
        for density_name in self.density_names:
            fit_out_subdir=os.path.join(fitout_dir, density_name+".mrc")
            prior_probability_file_path = os.path.join(fit_out_subdir, "prior_probabilities.txt")
            # 判断先验概率文件是否存在
            if os.path.exists(prior_probability_file_path):
                # 读取先验概率文件
                data=np.loadtxt(prior_probability_file_path, dtype=str)
                # 记录状态名
                states_of_densities.append([state for state in data[:,0]])
                states_to_id_dict.append({state:i for i, state in enumerate(data[:,0])})
                # 添加到先验概率列表
                prior_probs = [float(p) for p in data[:,1]]
                prior_prob_list.append(prior_probs)
                # 记录排序后的state_id
                prior_sorted_state_ids.append(np.argsort(prior_probs)[::-1].tolist())
            else:
                self.session.logger.error(f"Prior probability file {prior_probability_file_path} does not exist")
                return
        # 添加到self中
        self.states_of_densities = states_of_densities
        self.states_to_id_dict = states_to_id_dict
        self.prior_prob_list = prior_prob_list
        self.prior_sorted_state_ids = prior_sorted_state_ids
        # 标记先验结果已加载
        self.prior_results_loaded = True
        # 初始化状态选择器
        # 选择器记录先验结果中的state_id，用于更新显示
        prior_rank_list = [0 for i in range(len(self.density_names))]
        self.state_selection = [self.prior_sorted_state_ids[density_id][prior_rank_list[density_id]] for density_id in range(len(self.density_names))]

    # 获取后验概率结果
    def _get_posterior_results(self, project_directory, fitout_dir = ""):
        # 判断fitout_dir是否为空，为空则设置为默认值
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        # 验证fitout_dir是否存在
        if not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        # 后验概率列表
        posterior_prob_list = []
        # 记录每个状态的排名
        posterior_sorted_state_ids = []
        for density_id, density_name in enumerate(self.density_names):
            fit_out_subdir=os.path.join(fitout_dir, density_name+".mrc")
            posterior_probability_file_path = os.path.join(fit_out_subdir, "posterior_probabilities.txt")
            # 判断后验概率文件是否存在
            if os.path.exists(posterior_probability_file_path):
                # 读取后验概率文件
                data=np.loadtxt(posterior_probability_file_path, dtype=str)
                # 计算每一行的state_id，后续按次录入
                state_ids = [self.states_to_id_dict[density_id][state] for state in data[:,0]]
                index_list = np.argsort(state_ids)
                # 添加到后验概率列表
                posterior_probs = [float(p) for p in data[index_list,1]]
                posterior_prob_list.append(posterior_probs)
                # 记录按后验概率排序后的state_id
                posterior_sorted_state_ids.append(state_ids)
            else:
                self.session.logger.error(f"Posterior probability file {posterior_probability_file_path} does not exist")
                return
        # 添加到self中
        self.posterior_prob_list = posterior_prob_list
        self.posterior_sorted_state_ids = posterior_sorted_state_ids
        # 标记后验结果已加载
        self.posterior_results_loaded = True
        # 更新后验结果
        self._update_posterior_results()

        # 读取交联文件
        compliant_crosslinks_file_path = os.path.join(project_directory, "compliant_crosslinks.npy")
        if os.path.exists(compliant_crosslinks_file_path):
            self.compliant_crosslinks = np.load(compliant_crosslinks_file_path, allow_pickle=True).item()
        else:
            self.session.logger.error(f"Compliant crosslinks file {compliant_crosslinks_file_path} does not exist")
            return
        # 设置交联显示格式
        run(self.session,"distance style radius 0.3")
        # 如何没有compliant_crosslinks，则输出提示
        if len(self.compliant_crosslinks) == 0:
            self.session.logger.info("No compliant crosslinks found")
        else:
            # 绘制交联
            self._draw_crosslinks()
        

    # 更新当前状态到结果grid中
    def _update_current_states(self):
        for density_id, state_id in enumerate(self.state_selection):
            state = self.states_of_densities[density_id][state_id]
            doamin = "_".join(state.split("_")[:-1])
            fit_id = int(state.split("_")[-1])
            # 更新单密度状态
            # 更新单密度状态
            domain_name_text = self.result_layout.itemAtPosition(2+density_id, 1).widget()
            domain_name_text.setText(f"{doamin}")
            fit_id_text = self.result_layout.itemAtPosition(2+density_id, 2).widget()
            fit_id_text.setText(f"{fit_id}")
    
    # 导入并显示电镜密度
    def _open_density_files(self, map_directory):
        density_map_models=[]
        # 验证密度路径有效
        if not map_directory or not os.path.exists(map_directory) or not os.path.isdir(map_directory):
            self.session.logger.error(f"Map directory {map_directory} does not exist or is not a directory")
            return
        for density_name in self.density_names:
            density_file_path = os.path.join(map_directory, density_name+".mrc")
            if os.path.exists(density_file_path):
                map_model=run(self.session, f"open \"{density_file_path}\" name {density_name}")[0]
                map_model.set_parameters(surface_colors=[(178/255,178/255,178/255)],transparency=0.5)
                density_map_models.append(map_model)
            else:
                self.session.logger.error(f"Density file {density_file_path} does not exist")
                return
        # 调整视图
        run(self.session,"view")
        # 添加到self中
        self.density_map_models = density_map_models

    # 获取fit变换矩阵
    def get_transformation_matrix(self,log_path,fit_id):
        log_data=np.loadtxt(log_path,dtype=float,skiprows=fit_id,max_rows=1,usecols=range(2,14))
        transform_matrix=log_data.reshape((3,4))
        return transform_matrix
    
    # 模型平移旋转变换
    def transform_model(self,model,transform_matrix):
        model.position=Place(transform_matrix)
    
    # 导入当前状态的原子结构
    def _update_fitted_domains(self, project_directory, domain_directory = "", fitout_dir = ""):
        # 验证project_directory是否有效
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # 验证domain_directory是否为空，为空则设置为默认值
        if not domain_directory:
            domain_directory = os.path.join(project_directory, "domain_files")
        # 验证domain_directory是否有效
        if not os.path.exists(domain_directory):
            self.session.logger.error(f"Domain directory {domain_directory} does not exist")
            return
        # 验证fitout_dir是否为空，为空则设置为默认值
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        # 验证fitout_dir是否有效
        if not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        # 导入每个密度的状态
        for density_id, density_name in enumerate(self.density_names):
            fit_out_subdir=os.path.join(fitout_dir, density_name+".mrc")
            state_id = self.state_selection[density_id]
            state = self.states_of_densities[density_id][state_id]
            # 如果state发生变化，删除原有模型，导入新模型
            existing_state = self.fitted_domain_models[density_id][0]
            existing_model = self.fitted_domain_models[density_id][1]
            if state == existing_state:
                continue
            elif existing_state:
                # 关闭原有模型
                existing_model.delete()
            # 导入新模型
            domain = "_".join(state.split("_")[:-1])            
            fit_id = int(state.split("_")[-1])
            # 导入单密度状态
            domain_file_path = os.path.join(domain_directory, f"{domain}.pdb")
            if os.path.exists(domain_file_path):
                model = run(self.session, f"open \"{domain_file_path}\" name \"{density_name}: {domain}_{fit_id}\"")[0]
                # 导入fit log进行坐标变换
                log_path = os.path.join(fit_out_subdir, "fitlogs" , f"{domain}.log")
                transform_matrix = self.get_transformation_matrix(log_path,fit_id)
                self.transform_model(model,transform_matrix)
                # 更新记录
                self.fitted_domain_models[density_id] = (state, model)
            else:
                self.session.logger.error(f"Domain file {domain_file_path} does not exist")
                return
        
    # 根据当前状态在当前显示的模型中绘制交联
    def _draw_crosslinks(self, crosslink_color = "red", crosslink_radius = 1, crosslink_dashes = 5, label_height = 2.5):
        # 删除已有交联
        run(self.session,"distance delete")
        # 遍历所有交联
        count=0
        for density_1 in self.compliant_crosslinks.keys():
            density_id_1 = self.density_names.index(density_1)
            for state_1 in self.compliant_crosslinks[density_1].keys():
                current_state_id_1 = self.state_selection[density_id_1]
                current_state_1 = self.states_of_densities[density_id_1][current_state_id_1]
                if current_state_1 != state_1:
                    continue
                for density_2, state_2, residue_id_1, residue_id_2 in self.compliant_crosslinks[density_1][state_1]:
                    density_id_2 = self.density_names.index(density_2)
                    current_state_id_2 = self.state_selection[density_id_2]
                    current_state_2 = self.states_of_densities[density_id_2][current_state_id_2]
                    if current_state_2 != state_2:
                        continue
                    # 绘制交联
                    model_1 = self.fitted_domain_models[density_id_1][1]
                    model_id_1 = ".".join([str(item) for item in model_1.id])
                    model_2 = self.fitted_domain_models[density_id_2][1]
                    model_id_2 = ".".join([str(item) for item in model_2.id])
                    run(self.session,f"distance #{model_id_1}:{residue_id_1}@CA #{model_id_2}:{residue_id_2}@CA color {crosslink_color} radius {crosslink_radius} dashes {crosslink_dashes}")
                    count+=1
        # 设置标签格式
        if count > 0:
            run(self.session,f"label height {label_height}")
    


    # 更新先验结果到结果grid中
    def _update_prior_results(self):
        # 更新总体状态结果
        # aggregated_prior_prob_text = self.result_layout.itemAtPosition(1, 3).widget()
        # aggregated_prior_prob = np.prod([self.prior_prob_list[density_id][self.state_selection[density_id]] for density_id in range(len(self.density_names))])
        # aggregated_prior_prob_text.setText(f"{aggregated_prior_prob:.2e}")
        # 更新单密度状态结果
        for density_id, state_id in enumerate(self.state_selection):
            # 更新先验概率
            prior_prob_text = self.result_layout.itemAtPosition(2+density_id, 3).widget()
            prior_prob_text.setText(f"{self.prior_prob_list[density_id][state_id]:.4f}")
            # 更新单密度状态
            prior_rank = self.prior_sorted_state_ids[density_id].index(state_id)
            prior_rank_text = self.result_layout.itemAtPosition(2+density_id, 4).widget()
            prior_rank_text.setText(f"{prior_rank+1}")

    # 更新后验结果到结果grid中
    def _update_posterior_results(self):
        # 更新总体状态结果
        # 后验总概率不等于各概率之积，暂时不显示
        # 更新单密度状态结果
        for density_id, state_id in enumerate(self.state_selection):
            # 更新后验概率
            posterior_prob_text = self.result_layout.itemAtPosition(2+density_id, 5).widget()
            posterior_prob_text.setText(f"{self.posterior_prob_list[density_id][state_id]:.4f}")
            # 更新单密度状态
            posterior_rank = self.posterior_sorted_state_ids[density_id].index(state_id)
            posterior_rank_text = self.result_layout.itemAtPosition(2+density_id, 6).widget()
            posterior_rank_text.setText(f"{posterior_rank+1}")

    # 根据先验排名更新结果
    def _update_results_by_prior_ranks(self, project_directory, domain_directory = "", fitout_dir = ""):
        # 验证project_directory是否有效
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # 验证domain_directory是否为空，为空则设置为默认值
        if not domain_directory:
            domain_directory = os.path.join(project_directory, "domain_files")
        # 验证domain_directory是否有效
        if not os.path.exists(domain_directory):
            self.session.logger.error(f"Domain directory {domain_directory} does not exist")
            return
        # 验证fitout_dir是否为空，为空则设置为默认值
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        # 验证fitout_dir是否有效
        if not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        
        # 获取当前先验排名列表
        prior_rank_list = [int(prior_rank_text.text())-1 for prior_rank_text in [self.result_layout.itemAtPosition(2+density_id, 4).widget() for density_id in range(len(self.density_names))]]
        # 将排名转化为state_id
        state_id_list = [self.prior_sorted_state_ids[density_id][prior_rank] for density_id, prior_rank in enumerate(prior_rank_list)]
        # 更新当前状态选择器
        self.state_selection = state_id_list
        # 更新显示
        self._update_current_states()
        self._update_fitted_domains(project_directory, domain_directory, fitout_dir)
        if self.prior_results_loaded:
            self._update_prior_results()
        if self.posterior_results_loaded:
            self._update_posterior_results()
            # 绘制交联
            self._draw_crosslinks()

    # 根据后验排名更新结果
    def _update_results_by_posterior_ranks(self, project_directory, domain_directory = "", fitout_dir = ""):
        # 验证project_directory是否有效
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # 验证domain_directory是否为空，为空则设置为默认值
        if not domain_directory:
            domain_directory = os.path.join(project_directory, "domain_files")
        # 验证domain_directory是否有效
        if not os.path.exists(domain_directory):
            self.session.logger.error(f"Domain directory {domain_directory} does not exist")
            return
        # 验证fitout_dir是否为空，为空则设置为默认值
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        # 验证fitout_dir是否有效
        if not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        
        # 获取当前后验排名列表
        posterior_rank_list = [int(posterior_rank_text.text())-1 for posterior_rank_text in [self.result_layout.itemAtPosition(2+density_id, 6).widget() for density_id in range(len(self.density_names))]]
        # 将排名转化为state_id
        state_id_list = [self.posterior_sorted_state_ids[density_id][posterior_rank] for density_id, posterior_rank in enumerate(posterior_rank_list)]
        # 更新当前状态选择器
        self.state_selection = state_id_list
        # 更新显示
        self._update_current_states()
        self._update_fitted_domains(project_directory, domain_directory, fitout_dir)
        if self.prior_results_loaded:
            self._update_prior_results()
        if self.posterior_results_loaded:
            self._update_posterior_results()
            # 绘制交联
            self._draw_crosslinks()
                

    def _parse_domains(self, project_directory="",  pdb_dir="", pae_dir="", domain_directory="",\
                       plddt_cutoff=70, pae_cutoff=5, clique_cutoff=4, min_dege_ratio_between_cliques=0.6, min_common_nodes_ratio_between_cliques=0.5, minimum_domain_length=40, n_process=1):
        # Check if the project directory is valid, if not, raise an error
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # Check if the pdb directory is valid
        if not pdb_dir:
            pdb_dir = os.path.join(project_directory, "pdb_files").replace("\\", "/")
        elif not os.path.exists(pdb_dir):
            self.session.logger.error(f"PDB directory {pdb_dir} does not exist")
            return
        # Check if the pae directory is valid
        if not pae_dir:
            pae_dir = os.path.join(project_directory, "pae_files").replace("\\", "/")
        elif not os.path.exists(pae_dir):
            self.session.logger.error(f"PAE directory {pae_dir} does not exist")
            return
        # Check if the domain directory is valid
        if not domain_directory:
            output_dir = os.path.join(project_directory, "domain_files").replace("\\", "/")
        elif not os.path.exists(domain_directory):
            self.session.logger.error(f"Domain directory {domain_directory} does not exist")
            return
        else:
            output_dir = domain_directory.replace("\\", "/")
        # run parse_with_pae.py
        self.session.logger.info("Start parsing domains")
        # 打印参数
        self.session.logger.info(f"pdb directory: {pdb_dir}")
        self.session.logger.info(f"pae directory: {pae_dir}")
        self.session.logger.info(f"Output directory: {output_dir}")
        self.session.logger.info(f"PLDDT cutoff: {plddt_cutoff}")
        self.session.logger.info(f"PAE cutoff: {pae_cutoff}")
        self.session.logger.info(f"Clique cutoff: {clique_cutoff}")
        self.session.logger.info(f"Minimum dege ratio between cliques: {min_dege_ratio_between_cliques}")
        self.session.logger.info(f"Minimum common nodes ratio between cliques: {min_common_nodes_ratio_between_cliques}")
        self.session.logger.info(f"Minimum domain length: {minimum_domain_length}")
        self.session.logger.info(f"Number of process: {n_process}")
        # 启动进程
        arg_list = [f'{script_dir}/parse_with_pae.py',
                    self.error_log_path,
                    pdb_dir,
                    pae_dir,
                    output_dir,
                    n_process,
                    plddt_cutoff,
                    pae_cutoff,
                    clique_cutoff,
                    min_dege_ratio_between_cliques,
                    min_common_nodes_ratio_between_cliques,
                    minimum_domain_length]
        self.run_detatched_subprocess(arg_list)

    def _fit_and_score(self, project_directory, densities_dir, threshold, resolution, n_search, negtive_laplacian_cutoff, positive_laplacian_cutoff , n_process, domains_dir = "", fitout_dir = ""):
        # Check if the project directory is valid, if not, raise an error
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        
        # Check if the densities directory is valid, if not, raise an error
        if not densities_dir or not os.path.exists(densities_dir):
            self.session.logger.error(f"Densities directory {densities_dir} does not exist")
            return
        
        # Check if the domains directory is valid
        if not domains_dir:
            domains_dir = os.path.join(project_directory, "domain_files")
        elif not os.path.exists(domains_dir):
            self.session.logger.error(f"Domains directory {domains_dir} does not exist")
            return
        
        # Check if the fitout directory is valid
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        elif not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return

        # run fit_and_score.py
        print("Start fitting and scoring domains")
        print(f"Densities directory: {densities_dir}")
        print(f"Domains directory: {domains_dir}")
        print(f"Fitout directory: {fitout_dir}")
        print(f"Threshold: {threshold}")
        print(f"Resolution: {resolution}")
        print(f"Number of search: {n_search}")
        print(f"Negative laplacian cutoff: {negtive_laplacian_cutoff}")
        print(f"Positive laplacian cutoff: {positive_laplacian_cutoff}")
        print(f"Number of process: {n_process}")
        print()
        print() 
        # 启动进程
        arg_list=[f"{script_dir}/fit_with_chimerax.py",
                  self.error_log_path,
                  domains_dir,
                  densities_dir,
                  fitout_dir,
                  threshold,
                  resolution,
                  n_search,
                  negtive_laplacian_cutoff,
                  positive_laplacian_cutoff,
                  n_process]
        if sys.platform == 'win32':
            wait = False
        else:
            wait = True  # Linux: 需要设置wait = True，暂不确定什么原因
        self.run_detatched_subprocess(arg_list,wait=wait)
    
    def _calculate_prior_probability(self, project_directory, map_dir, fitout_dir, box_num, min_data_per_box, relative_density_cutoff, z_score_offset):
        # Check if the project directory is valid, if not, raise an error
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # Check if the map directory is valid, if not, raise an error
        if not map_dir or not os.path.exists(map_dir):
            self.session.logger.error(f"Map directory {map_dir} does not exist")
            return
        # Check if the fitout directory is valid
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        elif not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        # run calculate_prior_probability.py
        self.session.logger.info("Start calculating prior probabilities")
        arg_list = [f"{script_dir}/calculate_prior_probabilities.py",
                    self.error_log_path,
                    map_dir,
                    fitout_dir,
                    box_num,
                    min_data_per_box,
                    relative_density_cutoff,
                    z_score_offset]
        self.run_detatched_subprocess(arg_list)

    # 计算后验概率
    def _calculate_posterior_probability(self, project_directory,origin_domain_dir,map_dir,map_level,fitout_dir,acceptor_prior_probability_cutoff,donor_prior_probability_cutoff,evidence_strenth,crosslink_files):
        # Check if the project directory is valid, if not, raise an error
        if not project_directory or not os.path.exists(project_directory):
            self.session.logger.error(f"Project directory {project_directory} does not exist")
            return
        # Check if the origin domain directory is valid, if not, raise an error
        if not origin_domain_dir:
            origin_domain_dir = os.path.join(project_directory, "domain_files")
        elif not os.path.exists(origin_domain_dir):
            self.session.logger.error(f"Origin domain directory {origin_domain_dir} does not exist")
            return
        origin_domain_dir = origin_domain_dir.replace("\\", "/")
        # Check if the map directory is valid, if not, raise an error
        if not map_dir or not os.path.exists(map_dir):
            self.session.logger.error(f"Map directory {map_dir} does not exist")
            return
        # Check if the fitout directory is valid
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        elif not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return
        fitout_dir = fitout_dir.replace("\\", "/")
        # Check if the crosslink files are valid
        for crosslink_file in crosslink_files:
            if not os.path.exists(crosslink_file):
                self.session.logger.error(f"Crosslink file {crosslink_file} does not exist")
                return
        # run calculate_posterior_probability.py
        self.session.logger.info("Start calculating posterior probabilities")
        arg_list = [f"{script_dir}/calculate_posterior_probabilities.py",
                    self.error_log_path,
                    project_directory,
                    origin_domain_dir,
                    map_dir,
                    map_level,
                    fitout_dir,
                    acceptor_prior_probability_cutoff,
                    donor_prior_probability_cutoff,
                    evidence_strenth]
        arg_list+=crosslink_files
        self.run_detatched_subprocess(arg_list)


    def run_detatched_subprocess(self, arg_list, wait=False):
        """启动一个跨平台的后台进程"""
        kwargs = {
            'shell': True   # 目前tqdm的调用方式需要使用Shell
        }

        # get dir of executable chimerax
        chimerax_dir = os.path.dirname(os.path.realpath(sys.executable))

        # get current environment
        env = os.environ.copy()

        # add chimerax dir to PATH
        env['PATH'] = os.pathsep.join([env['PATH'], chimerax_dir])

        # 修改PYTHONPATH环境变量
        env['PYTHONPATH'] = os.pathsep.join(sys.path)

        # set executable python
        python_exe = [file_name for file_name in os.listdir(chimerax_dir) if file_name.startswith('python')][0]
        python_exe = os.path.join(chimerax_dir, python_exe)
        # python_exe = "python"

        # wait=True    # 调试时，等待进程结束，方便查看输出

        # 异步执行
        if not wait:
            if sys.platform == 'win32':
                kwargs['creationflags'] = subprocess.DETACHED_PROCESS   # 本来在GUI子系统下，设置shell=True时，应该不经过终端，直接启动shell。但这个参数设置后，会强制用终端启动，且会断开和主进程的管道。
            else:
                kwargs['start_new_session'] = True  # Linux: 分离进程组
            # 用列表传递python_exe和参数，避免路径中有空格出错
            proc = subprocess.Popen([python_exe] + arg_list, **kwargs, env=env)
        if wait:    # 后续考虑用其他方式输出报错信息，弹窗或输出到log
            if sys.platform == 'win32':
                pass
            else:
                kwargs['start_new_session'] = True  # Linux: 分离进程组
            kwargs['stdin'] = subprocess.PIPE
            kwargs['stdout'] = subprocess.PIPE
            kwargs['stderr'] = subprocess.PIPE
            # proc = subprocess.Popen(cmd, **kwargs, env=env)
            proc = subprocess.Popen([python_exe] + arg_list, **kwargs, env=env)
            # 获取输出
            stdout, stderr = proc.communicate()
            # 打印输出
            self.session.logger.info("Stdout:")
            self.session.logger.info(stdout.decode())
            self.session.logger.info("Stderr:")
            self.session.logger.info(stderr.decode())
            self.session.logger.info("Done.")