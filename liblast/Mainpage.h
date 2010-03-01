/** @mainpage LibFminer (libfminer)
 *
 * LibFminer
 *
 * This is the Fminer library, available at http://github.com/amaunz/libfminer/tree/master .<br>
 * The Fminer frontend application is available from http://github.com/amaunz/fminer/tree/master .<br>
 * You may download the scientific documentation from http://cs.maunz.de . The paper is entitled "Large Scale Graph Mining using Backbone Refinement Classes".<br>
 * Contact details are located at the end of this page.
 *
 *  @section sec0 News
 *
 * <i>29 Apr 2009</i>: The Backbone Refinement Class paper (co-authored by Christoph Helma and Stefan Kramer) has been accepted for the <b><a href="http://www.sigkdd.org/kdd/2009/papers.html">KDD 2009</a></b> conference on Data Mining and Knowledge Discovery (Jun 28 - Jul 1 2009 in Paris) for a presentation at the conference and inclusion in the conference proceedings.<p />
 * <i>30 Apr 2009</i>: The paper has been selected for oral presentation at <b><a href="http://www.cs.kuleuven.be/~dtai/ilp-mlg-srl/index.php?CONF=mlg&CONT=accepted">MLG 2009</a></b>.<p />
 * <i>10 Jul 2009</i>: <a href="http://doi.acm.org/10.1145/1557019.1557089" target="_blank">KDD conference proceedings are online</a>.
 *
 *  @section sec1 Abstract
 *
 *  We present a new approach to large-scale graph mining
 *  based on so-called backbone refinement classes. The method
 *  efficiently mines tree-shaped subgraph descriptors under minimum
 *  frequency and significance constraints, using classes
 *  of fragments to reduce feature set size and running times.
 *  The classes are defined in terms of fragments sharing a common
 *  backbone. The method is able to optimize structural
 *  inter-feature entropy as opposed to occurrences, which is
 *  characteristic for open or closed fragment mining. In the
 *  experiments, the proposed method reduces feature set sizes
 *  by >90 % and >30 % compared to complete tree mining and
 *  open tree mining, respectively. Evaluation using crossvalidation
 *  runs shows that their classification accuracy is similar
 *  to the complete set of trees but significantly better than
 *  that of open trees. Compared to open or closed fragment
 *  mining, a large part of the search space can be pruned due
 *  to an improved statistical constraint (dynamic upper bound
 *  adjustment), which is also confirmed in the experiments in
 *  lower running times compared to ordinary (static) upper
 *  bound pruning. Further analysis using large-scale datasets
 *  yields insight into important properties of the proposed descriptors,
 *  such as the dataset coverage and the class size
 *  represented by each descriptor. A final cross-validation run
 *  confirms that the novel descriptors render large training sets
 *  feasible which previously might have been intractable.
 *
 * \htmlonly
 *  <table align="center" border="0" width="610">
 *     <tr align="center">
 *         <td>
 *     <a href="http://www.maunz.de/salm.jpg">
 *     <img src="/salm.jpg" width="300" alt="Euclidean embedding of a Backbone Refinement Class Descriptors and Salmonella Mutagenicity data" title="Euclidean embedding of a Backbone Refinement Class Descriptors and Salmonella Mutagenicity data" />
 *      </a>
 *         </td>
 *     </tr>
 *     <tr>
 *         <td align="center">
 *     <small>Co-occurrence-based 2D embedding of molecules and backbone refinement class features showing close to perfect separation of target classes along top left to bottom right. (De)activating features are (red) green, (In)active instances (salmon) blue. Data: <a href="http://github.com/amaunz/cpdbdata/tree/master">CPDB salmonella mutagenicity</a>; Euclidean embedding: <a href="http://www.cs.kuleuven.be/~dtai/ilp-mlg-srl/index.php?CONF=ims&amp;CONT=wiki&amp;id=paper:ilp:33" target="_blank">Schulz <i>et. al,</i></a> </small>
 *         </td>
 *     </tr>
 *  </table>
 *  \endhtmlonly
 *
 *  \subsection ssec1 License
 *
 *   LibFminer is licensed under the terms of the GNU General Public License (GPL, see LICENSE). LibFminer is derived from (i.e., includes code from) the following project, licensed under GPL:
 * - Gaston: <a href="http://doi.acm.org/10.1145/1014052.1014134" target="_blank">Siegfried Nijssen and Joost Kok. A Quickstart in Frequent Structure Mining Can Make a Difference. Proceedings of the SIGKDD, 2004</a> (http://www.liacs.nl/home/snijssen/gaston/)
 *
 *   LibFminer uses (i.e., links to) the following projects, also licensed under GPL:
 * - <a href="http://openbabel.sourceforge.net/" target="_blank">OpenBabel</a>: The Open Babel Package, version 2.1.1.
 * - <a href="http://www.gnu.org/software/gsl/" target="_blank">GSL: GNU Scientific Library</a>, version 0.1.
 *
 *   These licensing conditions mean essentially, that your published program may only use (i.e., link to) and/or derive code from LibFminer under the condition that your source code is also freely available. Your personal usage of LibFminer, however, is not restricted in any way.
 *
 *  @section sec2 Using LibFminer
 *  LibFminer is a library, written in C++. It dynamically links to OpenBabel and GSL libraries.
 *  This section describes the installation of both the library and the frontend application for Linux and and 32 bit versions of Windows (NT and later).
 *  @subsection ssec20 Binary Quick Installation
 *  For Windows, you may use the <a href="http://github.com/amaunz/fminer/downloads" target="_blank">installer</a>. This installs the binary as well as libraries and C++ development headers. Otherwise install manually:
 *  - <a href="http://github.com/amaunz/libfminer/downloads" target="_blank">Download the binary DLL/SO file</a> and put it in a directory contained in your PATH/LD_LIBRARY_PATH environment variable.
 *  - <a href="http://github.com/amaunz/fminer/downloads" target="_blank">Download the binary frontend application</a>.
 *  - <b>Windows</b>: <a href="http://github.com/amaunz/openbabel-dll/tree" target="_blank">download the binary OpenBabel</a>. Click on the 'Download' button and chose 'ZIP'. After downloading, unpack and check file integrity (see README.txt). Then, put the DLL files in a directory contained in your PATH environment variable. <b>Linux</b>: install the package libopenbabel-dev.
 *  - <b>Windows</b>: <a href="http://github.com/amaunz/gsl-dll/tree" target="_blank">download the binary GSL</a>. Click on the 'Download' button and chose 'ZIP'. After downloading, unpack and check file integrity (see README.txt). Then, put the DLL files in a directory contained in your PATH environment variable. <b>Linux</b>: install the package libgsl0-dev.
 *
 *  @subsection ssec21 Compiling from source
 *  LibFminer is built as a dynamically loadable library.<br>
 *  Windows DLL: 
 *  - <a href="http://www.mingw.org/" target="_blank">Install Msys and MinGW</a>, then update gcc-core and gcc-g++ packages manually to the latest version.
 *  - OpenBabel: <a href="http://openbabel.org/wiki/Install_(MinGW)" target="_blank">follow the installation instrucations</a> to build OpenBabel yourself or <a href="http://github.com/amaunz/openbabel-dll/tree" target="_blank">download the binary DLLs</a> (with git: <code>git clone git://github.com/amaunz/openbabel-dll.git</code>).
 *  - GSL: <a href="http://github.com/amaunz/gsl-dll/tree" target="_blank">download the binary DLLs</a> (with git: <code>git clone git://github.com/amaunz/gsl-dll.git</code>).
 *  - <a href="http://github.com/amaunz/libfminer/tree" target="_blank">Download the library source code</a> (with git: <code>git clone git://github.com/amaunz/libfminer.git</code>). The <code>Makefile</code> automagically detects Windows. However, you have to adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - <a href="http://github.com/amaunz/fminer/tree" target="_blank">Download the frontend source code</a> (with git: <code>git clone git://github.com/amaunz/fminer.git</code>). In the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - To create this documentation with doxygen, type 'make doc'. The documentation explains API, constructor usage and options.
 *
 *  Linux SO: install development tools (gcc, g++, GNU make) and GSL as well as OpenBabel development package, then compile LibFminer. On Ubuntu, you can e.g. do it like this:
 *  - Install build tools and GSL:
 *    \code
 *    apt-get install build-essential             # development tools
 *    apt-get install libgsl0-dev                 # GSL binary lib and headers
 *    \endcode
 *  - OpenBabel: follow the installation instrucations at http://openbabel.org/wiki/Install_(source_code) to build yourself after doing:
 *    \code
 *    apt-get build-dep libopenbabel-dev          # build dependencies for OB
 *    apt-get source libopenbabel-dev             # extracts OB source to the current dir
 *    \endcode
 *    or try the repository version:
 *    \code
 *    apt-get install libopenbabel-dev            # OB binary lib and headers
 *    \endcode
 *  - <a href="http://github.com/amaunz/libfminer/tree" target="_blank">Download the library source code</a> (with git: <code>git clone git://github.com/amaunz/libfminer.git</code>). In the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - <a href="http://github.com/amaunz/fminer/tree" target="_blank">Download the frontend code</a> (with git: <code>git clone git://github.com/amaunz/fminer.git</code>). In the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - To create this documentation with doxygen, type 'make doc'. The documentation explains API, constructor usage and options.
 *  @subsection ssec22 Language Portability
 *  The API can be made available to other languages (currently only on Linux). Follow the installation instructions above. A config file for Swig to automagically create languages bindings exists (<code>fminer_wrap.i</code>). The Makefile also features a target that creates ruby bindings using this file. On Ubuntu, you can e.g. do this:
 *  - Swig: 
 *    \code
 *    apt-get install swig1.3 swig1.3-doc swig1.3-examples
 *    \endcode
 *  - Run <code>make ruby</code>.
 *  @section sec3 Pocket examples using the LibFminer API
 *  LibFminer uses the 'singleton' design pattern known from software engineering, i.e., class instantiation is restricted to one object. If you want to empty the database after a run to feed new compounds, use the Fminer::Reset() routine. 
 *  The following code demonstrate the use of the Fminer API from C++ and ruby. It feeds a set of class-labelled molecules in SMILES format (the API currently allows no gSpan input, use the frontend application for that) and calculates a vector of fragments along with statistical relevance and occurrences and prints them out. Every root node corresponds to a single chemical element. The output consists of gSpan graphs.
 *
 * Define the FMINER_SMARTS environment variable to produce output in SMARTS format. In this case, each line is a YAML sequence, containing SMARTS fragment, <i>p</i>-value, and two sequences denoting positive and negative class occurrences (line numbers in Smiles file): 
 *
 *  \code
 *  - [ smarts,    p_chisq,    occ_list_active,    occ_list_inactive ]
 *  \endcode
 *
 * Documentation for YAML can be found at: http://yaml.org/spec/cvs/current.html# Additionally define the FMINER_LAZAR environment variable to produce output in linfrag format which can be used as input to <code>Lazar</code>. 
 *
 *
 * \subsection CPP C++
 *
 * \code
 * #include "fminer.h"
 * #include <iostream>
 * #include <string.h>
 * using namespace std;
 *
 * Fminer* MyFminer;
 * int main(int argc, char *argv[], char *envp) {
 *   MyFminer= new Fminer();
 *   MyFminer->AddCompound ("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4" , 1);
 *   MyFminer->AddCompound ("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2" , 2);
 *      // ... continue adding compounds
 *   MyFminer->AddActivity((bool) true, 1);
 *   MyFminer->AddActivity((bool) false, 2);
 *      // ... continue adding activities (true for active, false for inactive)
 *   cerr << MyFminer->GetNoCompounds() << " compounds" << endl;
 *   // Toy example: special settings for mining all fragments
 *   MyFminer->SetChisqSig(0); // use no significance constraint
 *   MyFminer->SetRefineSingles(true); // refine structures with support 1
 *   // gather results for every root node in vector instead of immediate output
 *   MyFminer->SetConsoleOut(false);
 *   for ( int j = 0; j < (int) MyFminer->GetNoRootNodes(); j++ ) {
 *      vector<string>* result = MyFminer->MineRoot(j);
 *      for( int i = 0; i < result->size(); i++) {
 *        cout << (*result)[i] << endl;
 *      }
 *   }
 *   delete MyFminer;
 * }
 *
 *  \endcode
 *
 * \subsection Ruby Ruby
 *
 * This example assumes that you have created ruby bindings using <code>make fminer.so</code>.
 * \code
 *
 * require 'fminer'
 * MyFminer = Fminer::Fminer.new()
 * MyFminer.AddCompound("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4" , 1)
 * MyFminer.AddCompound("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2" , 2)
 *    # ... continue adding compounds
 * MyFminer.AddActivity(true, 1)
 * MyFminer.AddActivity(false, 2)
 *    # ... continue adding activities (true for active, false for inactive)
 * print MyFminer.GetNoCompounds()  
 * puts " compounds"
 * # Toy example: special settings for mining all fragments
 * # use no significance constraint
 * MyFminer.SetChisqSig(0) 
 * # refine structures with support 1
 * MyFminer.SetRefineSingles(true) 
 * # gather results for every root node in vector instead of immediate output
 * MyFminer.SetConsoleOut(false)
 * (0 .. MyFminer.GetNoRootNodes()-1).each do |j|
 *    result = MyFminer.MineRoot(j)
 *    puts "Results"
 *    (0 .. result.size-1).each do |i|
 *        puts result[i]
 *   end
 * end
 *
 *  \endcode
 *
 * \section Const Description of Constructors and Options
 * 
 * For the purpose of demonstration we used a toy database of two compounds and an unusual parameter configuration. Please note, that in general the defaults set by the standard constructor are sensible for most databases. They switch on BBRC mining with upper bound pruning for 95% significance and a minimum frequency of 2. The complete standard settings are:
 *
 * Minimum frequency: <b>2</b>, Feature type: <b>Trees</b>, Mine BBRCs: <b>true</b>, Dynamic upper bound: <b>true</b>, Significance level: <b>95%</b>, Console output: <b>false</b>, Aromatic perception: <b>false</b>, Refine Singles: <b>false</b>, Do Output: <b>true</b>, Separate BBRCs in output by blank line/vector: <b>false</b>, Most Specific Trees Only (Positive Border): <b>false</b>, Use line numbers instead of IDs: <b>false</b>
 *
 * \code
 *  //! Constructor for standard settings: 95% significance level, minimum frequency 2, type trees, dynamic upper bound, BBRC
 *  Fminer ();
 *  \endcode

 * There also exist more flexible constructors:
 * \code
 * //! Like standard constructor, but type and minimum frequency configurable
 * Fminer (int type, unsigned int minfreq);
 * //! Like standard constructor, but type, minimum frequency, significance level and BBRC configurable
 * Fminer (int type, unsigned int minfreq, float chisq_val, bool do_backbone);
 * \endcode
 * It is recommended to increase minimum frequency as a first step when too many features are generated. 

 * @section Contact Contact
 * Dipl.-Inf. Andreas Maunz<br>
 * Freiburg Center for Data Analysis and Modelling<br>
 * Hermann-Herder-Str. 3a<br>
 * 79104 Freiburg, Germany<br>
  Phone: +49761/203-8442, Fax: +49761/203-7700<br>
 * Email: maunza@fdm.uni-freiburg.de<br>
 * Web: http://cs.maunz.de
 *
 *  \author © 2008 by Andreas Maunz, 2008
 *
 * \htmlonly
 * <a href="http://www2.clustrmaps.com/counter/maps.php?url=http://www.maunz.de/libfminer-doc/main.html" id="clustrMapsLink"><img src="http://www2.clustrmaps.com/counter/index2.php?url=http://www.maunz.de/libfminer-doc/main.html" style="border:0px;visibility:hidden" alt="Locations of visitors to this page" title="Locations of visitors to this page" id="clustrMapsImg" />
 * </a>
 * <script type="text/javascript">
 * function cantload() {
 * img = document.getElementById("clustrMapsImg");
 * img.onerror = null;
 * img.src = "http://clustrmaps.com/images/clustrmaps-back-soon.jpg";
 * document.getElementById("clustrMapsLink").href = "http://clustrmaps.com";
 * }
 * img = document.getElementById("clustrMapsImg");
 * img.onerror = cantload;
 * </script>
 * \endhtmlonly
 **/
