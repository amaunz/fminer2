/** @mainpage LibBbrc (libbbrc)
 *
 * LibBbrc
 *
 * This is the Bbrc library, available at http://github.com/amaunz/fminer2/tree/master , subdirectory <code>libbbrc</code> (see below for download and build instructions).<br>
 * The Fminer frontend application is available from http://github.com/amaunz/fminer2/tree/master , subdirectory <code>fminer</code>.<br>
 * You may download the scientific documentation from http://cs.maunz.de . The paper is entitled "Large Scale Graph Mining using Backbone Refinement Classes".<br>
 * Supporting information is available here: http://bbrc.maunz.de .
 *
 *  @section Contents
 *  <ul>
 *      <li><a href="#News">News</a></li>
 *      <li><a href="#Abstract">Abstract</a></li>
 *      <li><a href="#Installation">Installation</a></li>
 *      <li><a href="#Guidance">Guidance on Using (Lib)Bbrc</a></li>
 *  </ul>
 * Contact details are located at the end of this page.
 *
 * <br><br>
 *  <a name="News">
 *  @section sec0 News
 *
 * <i>29 Apr 2009</i>: The Backbone Refinement Class paper (co-authored by Christoph Helma and Stefan Kramer) has been accepted for the <b><a href="http://www.sigkdd.org/kdd/2009/papers.html">KDD 2009</a></b> conference on Data Mining and Knowledge Discovery (Jun 28 - Jul 1 2009 in Paris) for a presentation at the conference and inclusion in the conference proceedings.<p />
 * <i>30 Apr 2009</i>: The paper has been selected for oral presentation at <b><a href="http://www.cs.kuleuven.be/~dtai/ilp-mlg-srl/index.php?CONF=mlg&CONT=accepted">MLG 2009</a></b>.<p />
 * <i>10 Jul 2009</i>: <a href="http://doi.acm.org/10.1145/1557019.1557089" target="_blank">KDD conference proceedings are online</a>.<p />
 * <i>29 Nov 2009</i>: Added configure script and bugfixes.<p />
 * <i>19 Apr 2010</i>: <a href="http://www.springerlink.com/content/c6m4w84794645213/" target="_blank">Machine Learning paper</a> accepted.<p />
 *
 * <br><br>
 *  <a name="Abstract">
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
 *     <small>Co-occurrence-based 2D embedding of molecules and backbone refinement class features showing close to perfect separation of target classes along top left to bottom right. (De)activating features are (red) green, (In)active instances (salmon) blue. Data: <a href="http://github.com/amaunz/cpdbdata/tree/master">CPDB salmonella mutagenicity</a>; Euclidean embedding: <a href="http://www.cs.kuleuven.be/~dtai/ilp-mlg-srl/index.php?CONF=ims&amp;CONT=wiki&amp;id=paper:ilp:33" target="_blank">Schulz <i>et. al,</i></a>.<br /> <b> <a href="/abstree.html" target="_blank" alt="Co-occurrence-based 2D embedding of molecules and backbone refinement class features">Click here</a> for a flash-animated version, indicating occurrences.</b></small>
 *         </td>
 *     </tr>
 *  </table>
 *  \endhtmlonly
 *
 *  \subsection ssec1 License
 *
 *   LibBbrc is licensed under the terms of the GNU General Public License (GPL, see LICENSE). LibBbrc is derived from (i.e. includes code from) the following project, licensed under GPL:
 * - <a href="http://doi.acm.org/10.1145/1014052.1014134" target="_blank">Siegfried Nijssen and Joost Kok. A Quickstart in Frequent Structure Mining Can Make a Difference. Proceedings of the SIGKDD, 2004</a> (http://www.liacs.nl/home/snijssen/gaston/)
 *
 *   LibBbrc uses (i.e. links to) the following projects, also licensed under GPL:
 * - <a href="http://openbabel.sourceforge.net/" target="_blank">OpenBabel</a>: The Open Babel Package, version 2.1.1.
 * - <a href="http://www.gnu.org/software/gsl/" target="_blank">GSL: GNU Scientific Library</a>, version 0.1.
 *
 *   These licensing conditions mean essentially that your published program may only use (i.e., link to) and/or derive code from LibBbrc under the condition that your source code is also freely available. This is to secure public availability and freedom of use.
 *
 * <br><br>
 *  <a name="Installation">
 *  @section sec2 Installation
 *  LibBbrc is a library, written in C++. It dynamically links to OpenBabel and GSL libraries.
 *  This section describes the installation of both the library and the frontend application for Linux.
 *  @subsection ssec22 Compiling from source
 *  <b>Linux SO</b>: install development tools (gcc, g++, GNU make) and GSL as well as OpenBabel development package, then compile LibBbrc. On Ubuntu, you can e.g. do it like this:
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
 *  - <a href="http://github.com/amaunz/fminer2/tree" target="_blank">Download the library source code</a> by clicking on "Download Source" (with git: <code>git clone git://github.com/amaunz/fminer2.git</code>) and cd to <code>libbbrc</code> subdirectory. Use <code>./configure</code> to configure the Makefile automatically or, in the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - Cd to <code>fminer</code> subdirectory. Use <code>./configure</code> to configure the Makefile automatically or, in the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - To create this documentation with doxygen, type 'make doc'. The documentation explains API, constructor usage and options.
 *  @subsection ssec23 Language Portability
 *  The API can be made available to other languages. Follow the installation instructions above. Config files for different languages bindings exist (<code>*bbrc_wrap.i</code>). 
 *
 *  The Makefile features a target that creates <b>ruby</b> bindings using this file. On Ubuntu, you can e.g. do this:
 *  - Use <code>./configure <version></code> to configure the Makefile automatically or, adjust the include flags (-I) in the Makefile in the line <code>INCLUDE_RB = ...</code> so that the directory contains file <code>ruby.h</code>. Also, let <code>RUBY = ...</code> point to the right ruby executable.
 *  - Run <code>make ruby</code>. Use <code>make rbtest</code> to test.
 *
 *  The Makefile features a target that creates <b>python</b> bindings using this file. On Ubuntu, you can e.g. do this:
 *  - Adjust the include flags (-I) in the Makefile in the line <code>INCLUDE_PY = ...</code> so that the directory contains file <code>Python.h</code>. Also, let <code>PYTHON = ...</code> point to the right ruby executable.
 *  - Run <code>make python</code>. Use <code>make pytest</code> to test.
 *
 *  <a name="Guidance">
 * @section Guidance Guidance on Using (Lib)Bbrc
 *
 * Most setting are sensible by default, see description of constructors and objects below. 
 *
 * I would suggest to manipulate the minimum frequency only at first. The number of fragments output should not be more than a multitude of the number of input graphs.
 * For minimum frequency, LibBbrc does not support percentage values. You will have to calculate absolute numbers.
 *
 *  @subsection sec33 Environment Variables
 *  <b>FMINER_SMARTS</b>: Produce output in SMARTS format. In this case, each line is a YAML sequence, containing SMARTS fragment, <i>p</i>-value, and two sequences denoting positive and negative class occurrences (line numbers in Smiles file): 
 *
 *  \code
 *  - [ smarts,    p_chisq,    occ_list_active,    occ_list_inactive ]
 *  \endcode
 *
 * Documentation for YAML can be found at: http://yaml.org/spec/cvs/current.html# (e.g. export FMINER_SMARTS=1).<br />
 * <b>FMINER_LAZAR</b>: Produce output in linfrag format which can be used as input to <code><a href="http://lazar.in-silico.de" target="_blank">Lazar</a></code> (e.g. export FMINER_LAZAR=1).<br />
 * <b>FMINER_PVALUES</b>: Produce p-values instead of chi-square values (e.g. export FMINER_PVALUES=1).<br />
 * <b>FMINER_NO_AROMATIC</b>: Disallow aromatic bonds on non-aromatic bonds (e.g. export FMINER_NO_AROMATIC=1).
 *
 *  @subsection sec3 Examples using the LibBbrc API
 *  LibBbrc uses the 'singleton' design pattern known from software engineering, i.e., class instantiation is restricted to one object. To empty the database after a run to feed new compounds, use the Bbrc::Reset() routine. 
 *
 *  The following code demonstrate the use of the Bbrc API from C++, python and ruby. It feeds a set of class-labelled molecules in SMILES format (the API currently allows no gSpan input, use the frontend application for that) and calculates a vector of fragments along with statistical relevance and occurrences and prints them out. Every root node corresponds to a single chemical element. The output consists of gSpan graphs.
 *
 *
 * \subsubsection CPP C++
 *
 * This example uses libBbrc in a C++ program.
 * The example assumes that you have created the C++ library using <code>make</code>.
 *
 * \code
 * #include "bbrc.h"
 * #include <iostream>
 * #include <string.h>
 * using namespace std;
 *
 * Bbrc* MyFminer;
 * int main(int argc, char *argv[], char *envp) {
 *   MyFminer= new Bbrc();
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
 *  \endcode
 *
 * \subsubsection Python Python
 *
 * This example assumes that you have created python bindings using <code>make python</code>.
 * \code
 * import bbrc
 * MyFminer = bbrc.Bbrc()
 * MyFminer.AddCompound("COC1=CC=C(C=C1)C2=NC(=C([NH]2)C3=CC=CC=C3)C4=CC=CC=C4" , 1)
 * MyFminer.AddCompound("O=C1NC(=S)NC(=O)C1C(=O)NC2=CC=CC=C2" , 2)
 * # ... continue adding compounds
 * MyFminer.AddActivity(1.0, 1)
 * MyFminer.AddActivity(0.0, 2)
 * # ... continue adding activities (true for active, false for inactive)
 * print repr(MyFminer.GetNoCompounds()) + ' compounds.'
 * # Toy example: special settings for mining all fragments
 * # use no significance constraint
 * MyFminer.SetChisqSig(0) 
 * # refine structures with support 1
 * MyFminer.SetRefineSingles(1) 
 * # gather results for every root node in vector instead of immediate output
 * MyFminer.SetConsoleOut(0)
 * for j in range(0, MyFminer.GetNoRootNodes()-1):
 *      result = MyFminer.MineRoot(j);
 *      for i in range(0, result.size()-1):
 *                  print result[i];
 * \endcode
 *
 * \subsubsection Ruby Ruby
 *
 * This example assumes that you have created ruby bindings using <code>make ruby</code>.
 * \code
 * require 'bbrc'
 * MyFminer = Bbrc::Bbrc.new()
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
 *    result.each do |res|
 *        puts res
 *   end
 * end
 *  \endcode
 *
 * \subsubsection Const Description of Constructors and Options
 * 
 * For the purpose of demonstration we used a toy database of two compounds and an unusual parameter configuration. Please note, that in general the defaults set by the standard constructor are sensible for most databases. They switch on BBRC mining with upper bound pruning for 95% significance and a minimum frequency of 2. The complete standard settings are:
 *
 * <ul>
 *  <li>Minimum frequency: <b>2</b></li>
 *  <li>Feature type: <b>Trees</b></li>
 *  <li>Mine BBRCs: <b>true</b></li>
 *  <li>Dynamic upper bound: <b>true</b></li>
 *  <li>Significance level: <b>95%</b></li>
 *  <li>Console output: <b>false</b></li>
 *  <li>Aromatic perception: <b>true</b></li>
 *  <li>Refine Singles: <b>false</b></li>
 *  <li>Do Output: <b>true</b></li>
 *  <li>Separate BBRCs in output by blank line/vector: <b>false</b></li>
 * </ul>
 *
 * \code
 *  //! Constructor for standard settings: 95% significance level, minimum frequency 2, type trees, dynamic upper bound, BBRC
 *  Bbrc ();
 *  \endcode

 * There also exist more flexible constructors:
 * \code
 * //! Like standard constructor, but type and minimum frequency configurable
 * Bbrc (int type, unsigned int minfreq);
 * //! Like standard constructor, but type, minimum frequency, significance level and BBRC configurable
 * Bbrc (int type, unsigned int minfreq, float chisq_val, bool do_backbone);
 * \endcode
 * It is recommended to increase minimum frequency as a first step when too many features are generated. 
 *
 * <br><br>
 * @section Contact Contact
 * Dipl.-Inf. Andreas Maunz<br>
 * Freiburg Center for Data Analysis and Modelling<br>
 * Hermann-Herder-Str. 3a<br>
 * 79104 Freiburg, Germany<br>
  Phone: +49761/203-8442, Fax: +49761/203-7700<br>
 * Email: maunza@fdm.uni-freiburg.de<br>
 * Web: http://cs.maunz.de
 *
 *  \author (c) 2010 by Andreas Maunz, 2010
 *
 **/
