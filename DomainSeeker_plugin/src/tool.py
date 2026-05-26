import os, sys
from chimerax.core.tools import ToolInstance
from chimerax.ui import MainToolWindow
from Qt.QtWidgets import (QVBoxLayout, QHBoxLayout, QGridLayout,
                          QWidget, QTabWidget, QSpacerItem,
                          QLabel,
                          QLineEdit, QTextEdit,
                          QPushButton, QCheckBox,
                          QFileDialog,
                          QFrame, QSizePolicy,
                          QScrollArea )
from Qt.QtCore import Qt
import numpy as np
from chimerax.core.commands import run
from chimerax.geometry.place import Place
import subprocess, threading

script_dir=os.path.dirname(os.path.realpath(__file__))

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
        calculation_tab = QWidget()
        calculation_tab_base_layout = QVBoxLayout(calculation_tab)

        scroll, calculation_layout = self._create_scroll_area("computation_scroll_area")
        calculation_tab_base_layout.addWidget(scroll)

        calculation_layout.addLayout(self._create_global_options_section())
        calculation_layout.addWidget(self.create_horizontal_line())
        calculation_layout.addLayout(self._create_file_fetching_section())
        calculation_layout.addWidget(self.create_horizontal_line())
        calculation_layout.addLayout(self._create_domain_parsing_section())
        calculation_layout.addWidget(self.create_horizontal_line())
        calculation_layout.addLayout(self._create_fitting_scoring_section())
        calculation_layout.addWidget(self.create_horizontal_line())
        calculation_layout.addLayout(self._create_prior_probability_section())
        calculation_layout.addWidget(self.create_horizontal_line())
        calculation_layout.addLayout(self._create_posterior_probability_section())
        calculation_layout.addWidget(self.create_horizontal_line(thickness=4))
        calculation_layout.addItem(QSpacerItem(0, 0, QSizePolicy.Minimum, QSizePolicy.Expanding))

        return calculation_tab

    def _create_global_options_section(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Global options"))
        body = QGridLayout()

        self.project_directory_text = QLineEdit()
        self._add_directory_row(body, 0, "Project directory:",
                                self.project_directory_text, read_only=True,
                                callback=lambda: self._initialize_project(self.project_directory_text))

        self.map_directory_text = QLineEdit()
        self._add_directory_row(body, 1, "Map directory:", self.map_directory_text,
                                callback=lambda: self._select_directory(
                                    self.map_directory_text, self.project_directory_text.text()))

        self.pdb_directory_text = QLineEdit()
        self._add_directory_row(body, 2, "Pdb directory:", self.pdb_directory_text,
                                callback=lambda: self._select_directory(
                                    self.pdb_directory_text, self.project_directory_text.text()))

        self.pae_directory_text = QLineEdit()
        self._add_directory_row(body, 3, "Pae directory:", self.pae_directory_text,
                                callback=lambda: self._select_directory(
                                    self.pae_directory_text, self.project_directory_text.text()))

        self.domain_directory_text = QLineEdit()
        self._add_directory_row(body, 4, "Domain directory:", self.domain_directory_text,
                                callback=lambda: self._select_directory(
                                    self.domain_directory_text, self.project_directory_text.text()))

        self.fitout_directory_text = QLineEdit()
        self._add_directory_row(body, 5, "Fitout directory:", self.fitout_directory_text,
                                callback=lambda: self._select_directory(
                                    self.fitout_directory_text, self.project_directory_text.text()))

        layout.addLayout(body)
        return layout

    def _create_file_fetching_section(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Fetch pdb and pae files from AFDB"))
        body = QHBoxLayout()
        body.addWidget(QLabel("Fetch proteins in file:"))

        protein_list_file_path_text = QLineEdit()
        body.addWidget(protein_list_file_path_text)

        select_btn = QPushButton("Select File")
        body.addWidget(select_btn)
        select_btn.clicked.connect(lambda: self._select_file(
            protein_list_file_path_text, self.project_directory_text.text()))

        fetch_btn = QPushButton("Fetch files")
        body.addWidget(fetch_btn)

        pae_only_cb = QCheckBox("PAE only")
        body.addWidget(pae_only_cb)

        fetch_btn.clicked.connect(lambda: self._fetch_pdb_and_pae_files(
            protein_list_file_path_text.text(), self.project_directory_text.text(),
            self.pdb_directory_text.text(), self.pae_directory_text.text(),
            pae_only_cb.isChecked()))
        layout.addLayout(body)
        return layout

    def _create_domain_parsing_section(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Parse proteins into domains based on PAE"))
        body = QGridLayout()

        plddt_cutoff_text = QLineEdit()
        plddt_cutoff_text.setText("70")
        body.addWidget(QLabel("plddt_cutoff"), 0, 0)
        body.addWidget(plddt_cutoff_text, 1, 0)

        pae_cutoff_text = QLineEdit()
        pae_cutoff_text.setText("5")
        body.addWidget(QLabel("pae_cutoff"), 0, 1)
        body.addWidget(pae_cutoff_text, 1, 1)

        clique_cutoff_text = QLineEdit()
        clique_cutoff_text.setText("4")
        body.addWidget(QLabel("clique_cutoff"), 0, 2)
        body.addWidget(clique_cutoff_text, 1, 2)

        min_dege_ratio_text = QLineEdit()
        min_dege_ratio_text.setText("0.6")
        body.addWidget(QLabel("min_dege_ratio"), 2, 0)
        body.addWidget(min_dege_ratio_text, 3, 0)

        min_common_nodes_ratio_text = QLineEdit()
        min_common_nodes_ratio_text.setText("0.5")
        body.addWidget(QLabel("min_common_nodes_ratio"), 2, 1, 1, 2)
        body.addWidget(min_common_nodes_ratio_text, 3, 1, 1, 2)

        min_domain_size_text = QLineEdit()
        min_domain_size_text.setText("40")
        body.addWidget(QLabel("min_domain_size"), 2, 3)
        body.addWidget(min_domain_size_text, 3, 3)

        max_domain_size_text = QLineEdit()
        max_domain_size_text.setText("1000")
        body.addWidget(QLabel("max_domain_size"), 2, 4)
        body.addWidget(max_domain_size_text, 3, 4)

        n_process_text = QLineEdit()
        n_process_text.setText("1")
        body.addWidget(QLabel("n_process"), 0, 3)
        body.addWidget(n_process_text, 1, 3)

        button = QPushButton("Parse domains")
        body.addWidget(button, 1, 4)
        button.clicked.connect(lambda: self._parse_domains(
            self.project_directory_text.text(), self.pdb_directory_text.text(),
            self.pae_directory_text.text(), self.domain_directory_text.text(),
            plddt_cutoff_text.text(), pae_cutoff_text.text(), clique_cutoff_text.text(),
            min_dege_ratio_text.text(), min_common_nodes_ratio_text.text(),
            min_domain_size_text.text(), max_domain_size_text.text(),
            n_process_text.text()))
        layout.addLayout(body)
        return layout

    def _create_fitting_scoring_section(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Fit domains into densities and score fitted domains"))
        body = QGridLayout()

        threshold_text = QLineEdit()
        threshold_text.setText("0.0")
        body.addWidget(QLabel("threshold"), 0, 0)
        body.addWidget(threshold_text, 1, 0)

        resolution_text = QLineEdit()
        resolution_text.setText("6.0")
        body.addWidget(QLabel("resolution"), 0, 1)
        body.addWidget(resolution_text, 1, 1)

        n_search_text = QLineEdit()
        n_search_text.setText("200")
        body.addWidget(QLabel("n_search"), 0, 2)
        body.addWidget(n_search_text, 1, 2)

        n_process_text = QLineEdit()
        n_process_text.setText("1")
        body.addWidget(QLabel("n_process"), 0, 3)
        body.addWidget(n_process_text, 1, 3)

        neg_cutoff_text = QLineEdit()
        neg_cutoff_text.setText("-0.001")
        body.addWidget(QLabel("negtive_laplacian_cutoff"), 2, 0)
        body.addWidget(neg_cutoff_text, 3, 0)

        pos_cutoff_text = QLineEdit()
        pos_cutoff_text.setText("0.001")
        body.addWidget(QLabel("positive_laplacian_cutoff"), 2, 1)
        body.addWidget(pos_cutoff_text, 3, 1)

        button = QPushButton("Fit & score")
        body.addWidget(button, 3, 3)
        button.clicked.connect(lambda: self._fit_and_score(
            self.project_directory_text.text(), self.map_directory_text.text(),
            threshold_text.text(), resolution_text.text(), n_search_text.text(),
            neg_cutoff_text.text(), pos_cutoff_text.text(), n_process_text.text(),
            self.domain_directory_text.text(), self.fitout_directory_text.text()))
        layout.addLayout(body)
        return layout

    def _create_prior_probability_section(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Calculate prior probability of each fitted domain"))
        body = QGridLayout()

        box_num_text = QLineEdit()
        box_num_text.setText("10")
        body.addWidget(QLabel("box_num"), 0, 0)
        body.addWidget(box_num_text, 1, 0)

        min_data_per_box_text = QLineEdit()
        min_data_per_box_text.setText("50")
        body.addWidget(QLabel("min_data_per_box"), 0, 1)
        body.addWidget(min_data_per_box_text, 1, 1)

        relative_density_cutoff_text = QLineEdit()
        relative_density_cutoff_text.setText("0.01")
        body.addWidget(QLabel("relative_density_cutoff"), 0, 2)
        body.addWidget(relative_density_cutoff_text, 1, 2)

        zScore_offset_text = QLineEdit()
        zScore_offset_text.setText("15")
        body.addWidget(QLabel("zScore_offset"), 2, 0)
        body.addWidget(zScore_offset_text, 3, 0)

        button = QPushButton("Calculate prior probability")
        body.addWidget(button, 3, 3)
        button.clicked.connect(lambda: self._calculate_prior_probability(
            self.project_directory_text.text(), self.map_directory_text.text(),
            self.fitout_directory_text.text(), box_num_text.text(),
            min_data_per_box_text.text(), relative_density_cutoff_text.text(),
            zScore_offset_text.text()))
        layout.addLayout(body)
        return layout

    def _create_posterior_probability_section(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Integrate extra experimental data"))
        body = QVBoxLayout()

        crosslink_layout = QVBoxLayout()
        crosslink_layout.addWidget(QLabel("XL-MS data"))

        sym_layout = QHBoxLayout()
        sym_layout.addWidget(QLabel("Surrounding symmetry config file:"))
        sym_text = QLineEdit()
        sym_layout.addWidget(sym_text)
        sym_btn = QPushButton("Select File")
        sym_layout.addWidget(sym_btn)
        sym_btn.clicked.connect(lambda: self._select_file(
            sym_text, self.project_directory_text.text()))
        crosslink_layout.addLayout(sym_layout)

        crosslink_body = QHBoxLayout()
        crosslink_files_manager = QVBoxLayout()

        select_btn = QPushButton("Select XL-MS files")
        crosslink_files_manager.addWidget(select_btn)

        crosslink_files_text = QTextEdit()
        crosslink_files_text.setReadOnly(True)
        crosslink_files_manager.addWidget(crosslink_files_text)
        select_btn.clicked.connect(lambda: self._select_files(
            crosslink_files_text, self.project_directory_text.text()))

        crosslink_body.addLayout(crosslink_files_manager, stretch=2)

        crosslink_options_layout = QGridLayout()

        post_threshold_text = QLineEdit()
        post_threshold_text.setText("0.0")
        crosslink_options_layout.addWidget(QLabel("threshold"), 0, 0)
        crosslink_options_layout.addWidget(post_threshold_text, 0, 1)

        acceptor_cutoff_text = QLineEdit()
        acceptor_cutoff_text.setText("0.00001")
        crosslink_options_layout.addWidget(QLabel("acceptor_cutoff"), 1, 0)
        crosslink_options_layout.addWidget(acceptor_cutoff_text, 1, 1)

        donor_cutoff_text = QLineEdit()
        donor_cutoff_text.setText("0.01")
        crosslink_options_layout.addWidget(QLabel("donor_cutoff"), 2, 0)
        crosslink_options_layout.addWidget(donor_cutoff_text, 2, 1)

        evidence_strength_text = QLineEdit()
        evidence_strength_text.setText("10")
        crosslink_options_layout.addWidget(QLabel("evidence_strength"), 3, 0)
        crosslink_options_layout.addWidget(evidence_strength_text, 3, 1)

        crosslink_body.addLayout(crosslink_options_layout, stretch=1)
        crosslink_layout.addLayout(crosslink_body)
        body.addLayout(crosslink_layout)

        button = QPushButton("Integrate experimental data")
        body.addWidget(button)
        button.clicked.connect(lambda: self._calculate_posterior_probability(
            self.project_directory_text.text(), self.domain_directory_text.text(),
            self.map_directory_text.text(), post_threshold_text.text(),
            self.fitout_directory_text.text(), acceptor_cutoff_text.text(),
            donor_cutoff_text.text(), evidence_strength_text.text(),
            sym_text.text(), crosslink_files_text.toPlainText().split("\n")))
        layout.addLayout(body)
        return layout
    
    # project 初始化
    def _initialize_project(self, project_directory_text):
        # 选择目录
        self._select_directory(project_directory_text)
        pass
        
                
        

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
        
        # 按钮：获取后验概率
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

    def _create_scroll_area(self, name, layout_class=QVBoxLayout):
        """Create a scroll area containing a widget with the given layout class."""
        scroll = QScrollArea()
        scroll.setObjectName(name)
        scroll.setWidgetResizable(True)
        container = QWidget()
        container.setObjectName(f"{name}_container")
        container.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        scroll.setWidget(container)
        layout = layout_class(container)
        return scroll, layout

    def _add_directory_row(self, layout, row, label, line_edit, read_only=False, callback=None):
        """Add a label + line edit + select button row to the given grid layout."""
        layout.addWidget(QLabel(label), row, 0)
        if read_only:
            line_edit.setReadOnly(True)
        layout.addWidget(line_edit, row, 1)
        button = QPushButton("Select Directory")
        layout.addWidget(button, row, 2)
        if callback:
            button.clicked.connect(callback)

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
                                                       "")
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
    
    def _fetch_pdb_and_pae_files(self, protein_list_file_path, project_directory, pdb_directory = "", pae_directory = "", pae_only = False):
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
        # Check if the pdb directory is valid (skip in pae_only mode)
        if not pae_only:
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
        self.session.logger.info("Start to fetch files...")
        if pae_only:
            arg_list=[f'{script_dir}/fetch_pdb_pae.py',
                      protein_list_file_path,
                      project_directory,
                      pae_directory,
                      '--pae-only']
        else:
            arg_list=[f'{script_dir}/fetch_pdb_pae.py',
                      protein_list_file_path,
                      project_directory,
                      pdb_directory,
                      pae_directory]
        self.run_detatched_subprocess(arg_list)

    def _reset_result_variables(self):
        """重置所有结果相关的运行时变量"""
        self.density_names = []
        self.states = []
        self.prior_probs = []
        self.posterior_idx = []
        self.posterior_probs = []
        self.state_selection = []
        self.fitted_domain_models = []
        self.density_map_models = {}
        self.compliant_crosslinks = {}
        self.symmetry_transform_list = []
        self.symmetry_models = {}
        self.prior_results_loaded = False
        self.posterior_results_loaded = False

    def _initialize_results(self, target_layout, map_directory, project_directory, domain_directory = "", fitout_dir = ""):
        # 重置运行时变量
        self._reset_result_variables()

        # 生成空白结果框架（复用已有布局则清空重填，否则新建）
        result_layout = self._generate_blank_results(target_layout, map_directory)
        self.result_layout = result_layout
        # 打开密度文件
        self._open_density_files(map_directory)
        # 初始化 fitted_domain_models 记录
        self.fitted_domain_models = [(None, None) for _ in range(len(self.density_names))]
        # 读取先验概率结果，并初始化选择器
        self._get_prior_results(project_directory, fitout_dir)
        # 显示初始状态
        self._update_current_states()
        # 显示初始状态的先验结果
        self._update_prior_results()
        # 更新原子模型
        self._update_fitted_domains(project_directory, domain_directory, fitout_dir)
        

        
    # 建立空白结果grid
    def _generate_blank_results(self, target_layout, map_directory):
        if hasattr(self, 'result_layout') and self.result_layout is not None:
            result_layout = self.result_layout
            # 清空布局中所有已有 widget
            while result_layout.count():
                item = result_layout.takeAt(0)
                widget = item.widget()
                if widget:
                    widget.deleteLater()
        else:
            result_scroll_area, result_layout = self._create_scroll_area("result_scroll_area", QGridLayout)
            result_layout.setObjectName("result_layout")
            result_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
            target_layout.addWidget(result_scroll_area)

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

        self.states = []
        self.prior_probs = []

        for density_name in self.density_names:
            fit_out_subdir=os.path.join(fitout_dir, density_name+".mrc")
            path = os.path.join(fit_out_subdir, "prior_probabilities.txt")
            if os.path.exists(path):
                data = np.loadtxt(path, dtype=str)
                self.states.append(data[:, 0].copy())
                self.prior_probs.append(data[:, 1].astype(np.float32))
            else:
                self.session.logger.error(f"Prior probability file {path} does not exist")
                return
        # 标记先验结果已加载
        self.prior_results_loaded = True
        # 初始化状态选择器（存 (state_name, prior_rank)，初始选最优先验排名）
        self.state_selection = [(self.states[d][0], 0) for d in range(len(self.density_names))]

    # 获取后验概率结果
    def _get_posterior_results(self, project_directory, fitout_dir = ""):
        # 判断fitout_dir是否为空，为空则设置为默认值
        if not fitout_dir:
            fitout_dir = os.path.join(project_directory, "fit_out")
        # 验证fitout_dir是否存在
        if not os.path.exists(fitout_dir):
            self.session.logger.error(f"Fitout directory {fitout_dir} does not exist")
            return

        self.posterior_idx = []
        self.posterior_probs = []

        for density_id, density_name in enumerate(self.density_names):
            fit_out_subdir=os.path.join(fitout_dir, density_name+".mrc")
            path = os.path.join(fit_out_subdir, "posterior_probabilities.txt")
            if os.path.exists(path):
                data = np.loadtxt(path, dtype=str)
                state_to_idx = {s: i for i, s in enumerate(self.states[density_id])}
                self.posterior_idx.append(
                    np.array([state_to_idx[s] for s in data[:, 0]], dtype=np.int32))
                self.posterior_probs.append(data[:, 1].astype(np.float32))
            else:
                self.session.logger.error(f"Posterior probability file {path} does not exist")
                return
        # 标记后验结果已加载
        self.posterior_results_loaded = True
        # 更新后验结果
        self._update_posterior_results()

        # 读取交联文件
        compliant_crosslinks_file_path = os.path.join(project_directory, "compliant_crosslinks.npy")
        if os.path.exists(compliant_crosslinks_file_path):
            compliant_crosslinks = np.load(compliant_crosslinks_file_path, allow_pickle=True).item()
        else:
            self.session.logger.error(f"Compliant crosslinks file {compliant_crosslinks_file_path} does not exist")
            return
        # 提取对称性变化列表到self.symmetry_transform_list
        # 如果compliant_crosslinks不存在"apply_symmetry_transform"键，设为空列表
        if "apply_symmetry_transform" not in compliant_crosslinks.keys():
            compliant_crosslinks["apply_symmetry_transform"] = []
        # 从self.compliant_crosslinks中提取对称性变化列表
        self.symmetry_transform_list = compliant_crosslinks["apply_symmetry_transform"]
        # 从compliant_crosslinks中删除"apply_symmetry_transform"键
        compliant_crosslinks.pop("apply_symmetry_transform")
        # 添加到self中
        self.compliant_crosslinks = compliant_crosslinks
        # 记录对称性模型{(density_name, copy_id):[density_model,domain_model]}
        self.symmetry_models = {}
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
        for density_id, (state_name, _) in enumerate(self.state_selection):
            domain = "_".join(state_name.split("_")[:-1])
            fit_id = int(state_name.split("_")[-1])
            domain_name_text = self.result_layout.itemAtPosition(2+density_id, 1).widget()
            domain_name_text.setText(f"{domain}")
            fit_id_text = self.result_layout.itemAtPosition(2+density_id, 2).widget()
            fit_id_text.setText(f"{fit_id}")
    
    # 导入并显示电镜密度
    def _open_density_files(self, map_directory):
        density_map_models={}
        # 验证密度路径有效
        if not map_directory or not os.path.exists(map_directory) or not os.path.isdir(map_directory):
            self.session.logger.error(f"Map directory {map_directory} does not exist or is not a directory")
            return
        for density_name in self.density_names:
            density_file_path = os.path.join(map_directory, density_name+".mrc")
            if os.path.exists(density_file_path):
                map_model=run(self.session, f"open \"{density_file_path}\" name {density_name}")[0]
                map_model.set_parameters(surface_colors=[(178/255,178/255,178/255)],transparency=0.5)
                density_map_models[density_name] = map_model
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
    
    # 从平移和旋转信息生成变换矩阵3*4
    def generate_transform_matrix(self, axis, center, angle_deg, translation):
        axis = np.array(axis, dtype=float)
        center = np.array(center, dtype=float)
        angle_deg = float(angle_deg)
        translation = np.array(translation, dtype=float)
        # 如果axis非空
        if len(axis)>0:
            # 计算旋转矩阵
            k = axis / np.linalg.norm(axis)
            angle_rad = np.deg2rad(angle_deg)
            cos_t, sin_t = np.cos(angle_rad), np.sin(angle_rad)

            K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
            R = np.eye(3) + sin_t * K + (1 - cos_t) * (K @ K)
        else:
            R = np.eye(3)
            center = np.zeros(3)
        # 如果translation是空
        if len(translation)==0:
            translation = np.zeros(3)
        t = (np.eye(3) - R) @ center + translation
        return np.hstack([R, t.reshape(3, 1)])
        
    
    # 模型变换
    def transform_model(self, model, transform_matrix):
        """Apply a 3x4 transform matrix [R|t] to the model's current position."""
        # 获取model自带的transform_matrix (3x4)
        M = model.position.matrix
        # 3x4 矩阵复合: new_tf @ M  =  [R_new @ R_M  |  R_new @ t_M + t_new]
        R_new, t_new = transform_matrix[:3, :3], transform_matrix[:3, 3]
        R_M, t_M = M[:3, :3], M[:3, 3]
        combined = np.hstack([R_new @ R_M, (R_new @ t_M + t_new).reshape(3, 1)])
        model.position = Place(combined)
    
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
            state, _ = self.state_selection[density_id]
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
        
    # 根据当前状态绘制交联
    def _draw_crosslinks(self, crosslink_color = "red", crosslink_radius = 1, crosslink_dashes = 5, label_height = 2.5):
        # 删除已有交联
        run(self.session,"distance delete")
        # 遍历所有交联
        count=0
        for density_1 in self.compliant_crosslinks.keys():
            density_id_1 = self.density_names.index(density_1)
            for state_1 in self.compliant_crosslinks[density_1].keys():
                current_state_1, _ = self.state_selection[density_id_1]
                if current_state_1 != state_1:
                    continue
                for item in self.compliant_crosslinks[density_1][state_1]:
                    density_2 = item[0]
                    state_2 = item[1]
                    residue_id_1 = item[2]
                    residue_id_2 = item[3]
                    density_id_2 = self.density_names.index(density_2)
                    current_state_2, _ = self.state_selection[density_id_2]
                    if current_state_2 != state_2:
                        continue
                    # 绘制交联
                    model_1 = self.fitted_domain_models[density_id_1][1]
                    model_id_1 = ".".join([str(item) for item in model_1.id])
                    copy_id_2 = 0 if len(item) < 5 else item[4]
                    # 如果copy_id_2>0，对应对称性模型，需要额外打开密度和原子模型，
                    # 并根据根据symmetry_transform_list[copy_id_2]进行对称性变换
                    if copy_id_2 > 0:
                        if (density_2, copy_id_2) not in self.symmetry_models.keys():
                            # init symmetry_model
                            self.symmetry_models[(density_2, copy_id_2)] = [None, None]
                            # symmetry matrix
                            symmetry_transform = self.symmetry_transform_list[copy_id_2-1]
                            axis = symmetry_transform["rotation_axis"]
                            center = symmetry_transform["rotation_point"]
                            angle_deg = symmetry_transform["rotation_degrees"]
                            translation = symmetry_transform["translation"]
                            symmetry_matrix = self.generate_transform_matrix(axis, center, angle_deg, translation)
                            # map
                            existing_map_model_2 = self.density_map_models[density_2]
                            # volume copy
                            symmetry_map_model_2 = run(self.session,f"volume copy #{existing_map_model_2.id[0]}")
                            # show orignal map
                            run(self.session,f"show #{existing_map_model_2.id[0]}")
                            # rename
                            run(self.session,f"rename #{symmetry_map_model_2.id[0]} {density_2}.{copy_id_2}")
                            # transform
                            self.transform_model(symmetry_map_model_2, symmetry_matrix)
                            # add to symmetry_models
                            self.symmetry_models[(density_2, copy_id_2)][0] = symmetry_map_model_2
                            # set color 浅黄色
                            symmetry_map_model_2.set_parameters(surface_colors=[(242/255,222/255,179/255)],transparency=0.5)
                            # domain
                            existing_domain_model_2 = self.fitted_domain_models[density_id_2][1]
                            # 复制：combine #{old_id} name {new_name}
                            # new_name格式："{density_name}.{copy_id}
                            symmetry_domain_model_2 = run(self.session,f" combine #{existing_domain_model_2.id[0]} name {density_2}.{copy_id_2}.pdb")
                            # transform
                            self.transform_model(symmetry_domain_model_2, symmetry_matrix)
                            # add to symmetry_models
                            self.symmetry_models[(density_2, copy_id_2)][1] = symmetry_domain_model_2
                            # 传递model
                            model_2 = symmetry_domain_model_2
                        else:
                            # 直接使用已有模型
                            symmetry_map_model_2, symmetry_domain_model_2 = self.symmetry_models[(density_2, copy_id_2)]
                            model_2 = symmetry_domain_model_2
                    else:
                        model_2 = self.fitted_domain_models[density_id_2][1]
                    model_id_2 = ".".join([str(item) for item in model_2.id])
                    run(self.session,f"distance #{model_id_1}:{residue_id_1}@CA #{model_id_2}:{residue_id_2}@CA color {crosslink_color} radius {crosslink_radius} dashes {crosslink_dashes}")
                    count+=1
        # 设置标签格式
        if count > 0:
            run(self.session,f"label height {label_height}")
    


    # 更新先验结果到结果grid中
    def _update_prior_results(self):
        # 更新单密度状态结果
        for density_id, (state_name, prior_rank) in enumerate(self.state_selection):
            prior_prob_text = self.result_layout.itemAtPosition(2+density_id, 3).widget()
            prior_prob_text.setText(f"{self.prior_probs[density_id][prior_rank]:.4f}")
            prior_rank_text = self.result_layout.itemAtPosition(2+density_id, 4).widget()
            prior_rank_text.setText(f"{prior_rank+1}")

    # 更新后验结果到结果grid中
    def _update_posterior_results(self):
        # 更新总体状态结果
        # 后验总概率不等于各概率之积，暂时不显示
        # 更新单密度状态结果
        for density_id, (state_name, prior_rank) in enumerate(self.state_selection):
            post_rank = int(np.where(self.posterior_idx[density_id] == prior_rank)[0][0])
            # 更新后验概率
            posterior_prob_text = self.result_layout.itemAtPosition(2+density_id, 5).widget()
            posterior_prob_text.setText(f"{self.posterior_probs[density_id][post_rank]:.4f}")
            # 更新单密度状态
            posterior_rank_text = self.result_layout.itemAtPosition(2+density_id, 6).widget()
            posterior_rank_text.setText(f"{post_rank+1}")

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
        # 排名 → (state 名, 先验排名)，更新状态选择器
        self.state_selection = [(self.states[d][prior_rank], prior_rank) for d, prior_rank in enumerate(prior_rank_list)]
        # 更新显示
        self._update_current_states()
        self._update_fitted_domains(project_directory, domain_directory, fitout_dir)
        if self.prior_results_loaded:
            self._update_prior_results()
        if self.posterior_results_loaded:
            self._update_posterior_results()
            # 删除所有对称性模型
            # 记录所有key，复制以防止迭代中修改字典导致错误
            keys = list(self.symmetry_models.keys()).copy()
            for key in keys:
                density_model, domain_model = self.symmetry_models[key]
                density_model.delete()
                domain_model.delete()
                # 从self.symmetry_models中移除该条目
                self.symmetry_models.pop(key)
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
        # 后验排名 → (state 名, 先验排名)，更新状态选择器
        self.state_selection = [(self.states[d][int(self.posterior_idx[d][posterior_rank])], int(self.posterior_idx[d][posterior_rank])) for d, posterior_rank in enumerate(posterior_rank_list)]
        # 更新显示
        self._update_current_states()
        self._update_fitted_domains(project_directory, domain_directory, fitout_dir)
        if self.prior_results_loaded:
            self._update_prior_results()
        if self.posterior_results_loaded:
            self._update_posterior_results()
            # 删除所有对称性模型
            # 记录所有key，复制以防止迭代中修改字典导致错误
            keys = list(self.symmetry_models.keys()).copy()
            for key in keys:
                density_model, domain_model = self.symmetry_models[key]
                density_model.delete()
                domain_model.delete()
                # 从self.symmetry_models中移除该条目
                self.symmetry_models.pop(key)
            # 绘制交联
            self._draw_crosslinks()
                

    def _parse_domains(self, project_directory="",  pdb_dir="", pae_dir="", domain_directory="",\
                       plddt_cutoff=70, pae_cutoff=5, clique_cutoff=4, \
                       min_dege_ratio_between_cliques=0.6, min_common_nodes_ratio_between_cliques=0.5, \
                       minimum_domain_length=40, maximum_domain_length=1000, \
                       n_process=1):
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
        # 启动进程
        self.session.logger.info("Start to parse domains...")
        arg_list = [f'{script_dir}/parse_with_pae.py',
                    pdb_dir,
                    pae_dir,
                    output_dir,
                    n_process,
                    plddt_cutoff,
                    pae_cutoff,
                    clique_cutoff,
                    min_dege_ratio_between_cliques,
                    min_common_nodes_ratio_between_cliques,
                    minimum_domain_length,
                    maximum_domain_length]
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
        self.session.logger.info("Start to fit and score domains...")
        # 启动进程
        arg_list=[f"{script_dir}/fit_with_chimerax.py",
                  domains_dir,
                  densities_dir,
                  fitout_dir,
                  threshold,
                  resolution,
                  n_search,
                  negtive_laplacian_cutoff,
                  positive_laplacian_cutoff,
                  n_process]
        self.run_detatched_subprocess(arg_list)
    
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
        self.session.logger.info("Start to calculate prior probabilities...")
        arg_list = [f"{script_dir}/calculate_prior_probabilities.py",
                    map_dir,
                    fitout_dir,
                    box_num,
                    min_data_per_box,
                    relative_density_cutoff,
                    z_score_offset]
        self.run_detatched_subprocess(arg_list)

    # 计算后验概率
    def _calculate_posterior_probability(self, project_directory,origin_domain_dir,map_dir,map_level,fitout_dir,acceptor_prior_probability_cutoff,donor_prior_probability_cutoff,evidence_strenth,symmetry_transform_file,crosslink_files):
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
        # Check if the symmetry transform file is valid
        # 如果为空，设置为“None”
        if not symmetry_transform_file:
            symmetry_transform_file = "None"
        elif not os.path.exists(symmetry_transform_file):
            self.session.logger.error(f"Symmetry transform file {symmetry_transform_file} does not exist")
            return
        # Check if the crosslink files are valid
        for crosslink_file in crosslink_files:
            if not os.path.exists(crosslink_file):
                self.session.logger.error(f"Crosslink file {crosslink_file} does not exist")
                return
        # run calculate_posterior_probability.py
        self.session.logger.info("Start to calculate posterior probabilities...")
        arg_list = [f"{script_dir}/calculate_posterior_probabilities.py",
                    project_directory,
                    origin_domain_dir,
                    map_dir,
                    map_level,
                    fitout_dir,
                    acceptor_prior_probability_cutoff,
                    donor_prior_probability_cutoff,
                    evidence_strenth,
                    symmetry_transform_file]
        arg_list+=crosslink_files
        self.run_detatched_subprocess(arg_list)

    # ── 子进程消息路由 ──────────────────────────────────────────────
    # stdout → info / status（[STATUS] 前缀 → 状态栏，其余 → info）
    # stderr → warning / error（[WARN] 前缀 → 警告，其余含 traceback → error）
    # ChimeraX 只允许主线程更新 UI，故 logger 调用经 thread_safe 调度

    def _read_stdout(self, stream):
        """读取子进程 stdout"""
        try:
            for line in iter(stream.readline, ''):
                line = line.rstrip()
                if line.startswith("[STATUS]"):
                    self.session.ui.thread_safe(self.session.logger.status, line[8:])
                elif line.startswith("[INFO]"):
                    self.session.ui.thread_safe(self.session.logger.info, line[6:])
                elif line:
                    self.session.ui.thread_safe(self.session.logger.info, line)
        except (ValueError, RuntimeError):
            pass

    def _read_stderr(self, stream):
        """读取子进程 stderr"""
        try:
            for line in iter(stream.readline, ''):
                line = line.rstrip()
                if line.startswith("[WARN]"):
                    self.session.ui.thread_safe(self.session.logger.warning, line[6:])
                elif line.startswith("[ERROR]"):
                    self.session.ui.thread_safe(self.session.logger.error, line[7:])
                else:
                    self.session.ui.thread_safe(self.session.logger.error, line)
        except (ValueError, RuntimeError):
            pass


    def run_detatched_subprocess(self, arg_list):
        """启动一个跨平台的后台进程"""
        # get dir of executable chimerax
        chimerax_dir = os.path.dirname(os.path.realpath(sys.executable))
        # mac特殊处理
        if sys.platform == 'darwin':
            chimerax_dir = chimerax_dir.replace("/Contents/MacOS", "/Contents/bin")

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

        kwargs = {
            'shell': False,
            "text": True,    # 输出为文本
            "bufsize": 1,   # 行缓冲
            "stdout": subprocess.PIPE,
            "stderr": subprocess.PIPE,
        }

        # 设置各系统参数
        # windows
        if sys.platform in ["win32","win64"]:
            # kwargs['creationflags'] = subprocess.DETACHED_PROCESS   # 本来在GUI子系统下，设置shell=True时，应该不经过终端，直接启动shell。但这个参数设置后，会强制用终端启动，且会断开和主进程的管道。
            kwargs['creationflags'] = subprocess.CREATE_NO_WINDOW   # Windows: 不显示窗口。后续捕获输出到chimerax的logger中
        # mac
        elif sys.platform == "darwin":
            pass
        # linux
        else:
            # kwargs['start_new_session'] = True  # Linux: 分离进程组
            # linux系统下，shell=True时，如果传递命令列表，会无参数启动第一个命令。因此，linux系统下，如果指定shell=True时，传递命令字符串
            pass

        # 启动进程
        proc = subprocess.Popen([python_exe] + arg_list, **kwargs, env=env) # 用列表传递python_exe和参数，避免路径中有空格出错
        # 监控 stdout（info / status）
        stdout_thread = threading.Thread(target=self._read_stdout, args=(proc.stdout,))
        stdout_thread.daemon = True
        stdout_thread.start()
        # 监控 stderr（warning / error）
        stderr_thread = threading.Thread(target=self._read_stderr, args=(proc.stderr,))
        stderr_thread.daemon = True
        stderr_thread.start()