/** @mainpage LibLast (liblast)
 *
 * LibLast
 *
 * This is the Last library, available at http://github.com/amaunz/fminer2/tree/master , subdirectory <code>liblast</code> (see below for download and build instructions).<br>
 * The Fminer frontend application is available from http://github.com/amaunz/fminer2/tree/master , subdirectory <code>fminer</code>.<br>
 * Supporting information is available here: http://last-pm.maunz.de .
 *
 *  @section Contents
 *  <ul>
 *      <li><a href="#News">News</a></li>
 *      <li><a href="#Abstract">Abstract</a></li>
 *      <li><a href="#Installation">Installation</a></li>
 *      <li><a href="#Guidance">Guidance on Using (Lib)Last</a></li>
 *  </ul>
 * Contact details are located at the end of this page.
 *
 * <br><br>
 *  <a name="News">
 *  @section sec0 News
 * <i>22 Jun 2010</i>: The LAST-PM paper (co-authored by Christoph Helma, Tobias Cramer and Stefan Kramer) has been accepted for the <b><a href="http://www.ecmlpkdd2010.org" target="_blank">ECML/PKDD 2010</a></b> conference (Sep 20 - Sep 24 2010 in Barcelona) for a presentation at the conference and inclusion in the conference proceedings.<p />
 *
 * <br><br>
 *  <a name="Abstract">
 *  @section sec1 Abstract
 *
 * Pattern mining methods for graph data have largely been restricted to ground features, such as frequent or correlated subgraphs. Kazius et al. have demonstrated the use of elaborate patterns in the biochemical domain, summarizing several ground features at once. Such patterns bear the potential to reveal latent information not present in any individual ground feature. However, those patterns were handcrafted by chemical experts. 
 * In this paper, we present a data-driven bottom-up method for pattern generation that takes advantage of the embedding relationships among individual ground features. The method works fully automatically and does not require data preprocessing (e.g., to introduce abstract node or edge labels). Controlling the process of generating ground features, it is possible to align them canonically and merge (stack) them, yielding a weighted edge graph. In a subsequent step, the subgraph features are compressed by singular value decomposition (SVD). Our experiments show that the resulting features are chemically meaningful and that they can enable substantial performance improvements on chemical datasets that have been problematic so far for graph mining approaches.
 *
 *  \subsection ssec1 License
 *
 *   LibLast is licensed under the terms of the GNU General Public License (GPL, see LICENSE). LibLast is derived from (i.e. includes code from) the following project, licensed under GPL:
 * - <a href="http://doi.acm.org/10.1145/1014052.1014134" target="_blank">Siegfried Nijssen and Joost Kok. A Quickstart in Frequent Structure Mining Can Make a Difference. Proceedings of the SIGKDD, 2004</a> (http://www.liacs.nl/home/snijssen/gaston/)
 *
 *   LibLast uses (i.e. links to) the following projects, also licensed under GPL:
 * - <a href="http://openbabel.sourceforge.net/" target="_blank">OpenBabel</a>: The Open Babel Package, version 2.1.1.
 * - <a href="http://www.gnu.org/software/gsl/" target="_blank">GSL: GNU Scientific Library</a>, version 0.1.
 *
 *   These licensing conditions mean essentially that your published program may only use (i.e., link to) and/or derive code from LibLast under the condition that your source code is also freely available. This is to secure public availability and freedom of use.
 *
 * <br><br>
 *  <a name="Installation">
 *  @section sec2 Installation
 *  LibLast is a library, written in C++. It dynamically links to OpenBabel and GSL libraries.
 *  This section describes the installation of both the library and the frontend application for Linux.
 *  @subsection ssec22 Compiling from source
 *  <b>Linux SO</b>: install development tools (gcc, g++, GNU make) and GSL as well as OpenBabel development package, then compile LibLast. On Ubuntu, you can e.g. do it like this:
 *  - Install build tools and GSL:
 *    \code
 *    apt-get install build-essential             # development tools
 *    apt-get install libgsl0-dev                 # GSL binary lib and headers
 *    \endcode
 *  - OpenBabel: follow the <a href="http://openbabel.org/wiki/Get_Open_Babel">installation instrucations</a> to build yourself after doing:
 *    \code
 *    apt-get build-dep libopenbabel-dev          # build dependencies for OB
 *    apt-get source libopenbabel-dev             # extracts OB source to the current dir
 *    \endcode
 *    or try the repository version:
 *    \code
 *    apt-get install libopenbabel-dev            # OB binary lib and headers
 *    \endcode
 *    Note: you will need the former if you want to use <a href="https://github.com/amaunz/last-utils">LAST-UTILS</a>, see the <a href="https://github.com/amaunz/last-utils/blob/master/INSTALL">INSTALL there</a>.
 *  - <a href="http://github.com/amaunz/fminer2/tree" target="_blank">Download the library source code</a> by clicking on "Download Source". (with git: <code>git clone git://github.com/amaunz/fminer2.git</code>) and cd to <code>liblast</code> subdirectory. Use <code>./configure</code> to configure the Makefile automatically or, in the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - Cd to <code>fminer</code> subdirectory. Use <code>./configure</code> to configure the Makefile automatically or, in the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - To create this documentation with doxygen, type 'make doc'. The documentation explains API, constructor usage and options.
 *  @subsection ssec23 Language Portability
 *  The API can be made available to other languages.
 *
 *  The Makefile features a target that creates <b>ruby</b> bindings. On Ubuntu, you can e.g. do this:
 *  - Use <code>./configure <version></code> to configure the Makefile automatically or, adjust the include flags (-I) in the Makefile in the line <code>INCLUDE_RB = ...</code> so that the directory contains file <code>ruby.h</code>. Also, let <code>RUBY = ...</code> point to the right executable.
 *  - Run <code>make ruby</code>. Use <code>make rbtest</code> to test. The configuration was tested with swig 1.3.40 and Ruby 1.8.
 *
 *  The Makefile features a target that creates <b>python</b> bindings. On Ubuntu, you can e.g. do this:
 *  - Adjust the include flags (-I) in the Makefile in the line <code>INCLUDE_PY = ...</code> so that the directory contains file <code>Python.h</code>. Also, let <code>PYTHON = ...</code> point to the right executable.
 *  - Run <code>make python</code>. Use <code>make pytest</code> to test. The configuration was tested with swig 1.3.40 and Python 2.5.
 *
 * <b>Important:</b> There are swig interface files (<code>*.i</code>) and pre-configured swig output files (<code>*.cxx</code>). You need to re-create those output files if you are deploying for newer versions of the target languages, and you can find the necessary swig calls in the Makefile (commented out).
 *
 *  <a name="Guidance">
 * @section Guidance Guidance on Using (Lib)Last
 *
 * Most setting are sensible by default, see description of constructors and objects below. 
 *
 * I would suggest to manipulate the minimum frequency only at first. The number of fragments output should not be more than a multitude of the number of input graphs.
 * For minimum frequency, LibLast does not support percentage values. You will have to calculate absolute numbers.
 *
 *  @subsection sec3 Examples using the LibLast API
 *  LibLast uses the 'singleton' design pattern known from software engineering, i.e., class instantiation is restricted to one object. To empty the database after a run to feed new compounds, use the Last::Reset() routine. 
 *
 *  The following code demonstrate the use of the Last API from C++, python, and ruby. It feeds a set of class-labelled molecules in SMILES format (the API currently allows no gSpan input, use the frontend application for that) and calculates a set of latent fragments and prints them out. Every root node corresponds to a single chemical element. The output consists of <a href="http://graphml.graphdrawing.org">GraphML</a> which can be postprocessed to SMARTS patterns using the <a href="http://github.com/amaunz/last-utils" target="_blank">LAST-UTILS</a>.
 *
 *  @subsection sec33 Environment Variables
 * <b>FMINER_SILENT</b>      : Redirect STDERR (debug output) of fminer to local file 'fminer_debug.txt'<br />
 *
 * Note: The value you set the environment variables to is irrelevant. Use <code>unset</code> to disable the environment variables, e.g. <code>unset FMINER_LAZAR</code>.
 *
 * \subsubsection CPP C++
 *
 * This example uses libLast in a C++ program.
 * The example assumes that you have created the C++ library using <code>make</code>.
 *
 * \code
 *
 * #include "last.h"
 * #include <iostream>
 * #include <string.h>
 * using namespace std;
 *
 * Last* MyFminer;
 * int main(int argc, char *argv[], char *envp[]) {
 *     MyFminer= new Last();
 *     MyFminer = Last::Last.new();
 *     MyFminer->SetMaxHops(25);
 *     MyFminer->SetConsoleOut(false);
 *     // Add compounds below. IMPORTANT! Do not change settings after adding compounds!
 *     MyFminer->AddCompound ("O=C(C(C(C=C3)=CC=C3O)=CO2)C1=C2C=C(O)C=C1O" , 1);
 *     MyFminer->AddCompound ("Oc1ccc(cc1)[C@@H]2Cc3ccc(O)cc3OC2" , 2);
 *     MyFminer->AddCompound ("O=C1C(C3=CC=C(O)C=C3)=COC2=C1C=CC(O)=C2" , 3);
 *     MyFminer->AddCompound ("O=C1C(C3=CC=C(OC)C=C3)=COC2=C1C=CC(O)=C2" , 4);
 *     MyFminer->AddCompound ("OC1=CC=C(CCCCCCCC)C=C1" , 5);
 *     MyFminer->AddCompound ("C1(C=CC=CC=1C(=C(Cl)Cl)C2=CC=C(C=C2)Cl)Cl" , 6);
 *     MyFminer->AddCompound ("O=C(C1=C(C=CC=C1)C(=O)OCC(CCCC)CC)OCC(CCCC)CC" , 7);
 *     MyFminer->AddCompound ("Oc1cc(O)cc2CCCCC[C@@H](O)CCC[C@H](C)OC(=O)c12" , 8);
 *     MyFminer->AddCompound ("O=C1C2=C(C=C(C=C2O)O)OC(=C1O)C3=CC(=C(C=C3)O)O" , 9);
 *     MyFminer->AddCompound ("C1(=C(C(=O)C2=C(O1)C=C(C=C2)O)O)C3=CC(O)=C(C=C3)O" , 10);
 *     //... continue adding compounds
 *     MyFminer->AddActivity((bool) 1.0, 1);
 *     MyFminer->AddActivity((bool) 1.0, 2);
 *     MyFminer->AddActivity((bool) 1.0, 3);
 *     MyFminer->AddActivity((bool) 1.0, 4);
 *     MyFminer->AddActivity((bool) 1.0, 5);
 *     MyFminer->AddActivity((bool) 1.0, 6);
 *     MyFminer->AddActivity((bool) 1.0, 7);
 *     MyFminer->AddActivity((bool) 1.0, 8);
 *     MyFminer->AddActivity((bool) 0.0, 9);
 *     MyFminer->AddActivity((bool) 0.0, 10);
 *     //... continue adding activities (1.0 for active, 0.0 for inactive)
 *     cerr << MyFminer->GetNoCompounds() << " compounds" << endl;
 *     // gather results for every root node in vector instead of immediate output
 *     for ( int j = 0; j < (int) MyFminer->GetNoRootNodes(); j++ ) {
 *        vector<string>* result = MyFminer->MineRoot(j);
 *        for ( int i = 0; i < result->size(); i++) {
 *           cout << (*result)[i] << endl;
 *        }
 *     }
 *     delete MyFminer;
 * }
 *
 * 
 * \endcode
 *
 * \subsubsection Python Python
 *
 * This example assumes that you have created python bindings using <code>make python</code>.
 *
 * \code
 * import liblast
 * MyFminer = liblast.Last()
 * MyFminer.SetMaxHops(25)
 * MyFminer.SetConsoleOut(0)
 * # Add compounds below. IMPORTANT! Do not change settings after adding compounds!
 * MyFminer.AddCompound("O=C(C(C(C=C3)=CC=C3O)=CO2)C1=C2C=C(O)C=C1O" , 1)
 * MyFminer.AddCompound("Oc1ccc(cc1)[C@@H]2Cc3ccc(O)cc3OC2" , 2)
 * MyFminer.AddCompound("O=C1C(C3=CC=C(O)C=C3)=COC2=C1C=CC(O)=C2" , 3)
 * MyFminer.AddCompound("O=C1C(C3=CC=C(OC)C=C3)=COC2=C1C=CC(O)=C2" , 4)
 * MyFminer.AddCompound("OC1=CC=C(CCCCCCCC)C=C1" , 5)
 * MyFminer.AddCompound("C1(C=CC=CC=1C(=C(Cl)Cl)C2=CC=C(C=C2)Cl)Cl" , 6)
 * MyFminer.AddCompound("O=C(C1=C(C=CC=C1)C(=O)OCC(CCCC)CC)OCC(CCCC)CC" , 7)
 * MyFminer.AddCompound("Oc1cc(O)cc2CCCCC[C@@H](O)CCC[C@H](C)OC(=O)c12" , 8)
 * MyFminer.AddCompound("O=C1C2=C(C=C(C=C2O)O)OC(=C1O)C3=CC(=C(C=C3)O)O" , 9)
 * MyFminer.AddCompound("C1(=C(C(=O)C2=C(O1)C=C(C=C2)O)O)C3=CC(O)=C(C=C3)O" , 10)
 * # ... continue adding compounds
 * MyFminer.AddActivity(1.0, 1)
 * MyFminer.AddActivity(1.0, 2)
 * MyFminer.AddActivity(1.0, 3)
 * MyFminer.AddActivity(1.0, 4)
 * MyFminer.AddActivity(1.0, 5)
 * MyFminer.AddActivity(1.0, 6)
 * MyFminer.AddActivity(1.0, 7)
 * MyFminer.AddActivity(1.0, 8)
 * MyFminer.AddActivity(0.0, 9)
 * MyFminer.AddActivity(0.0, 10)
 * # ... continue adding activities (true for active, false for inactive)
 * print repr(MyFminer.GetNoCompounds()) + ' compounds'
 * # gather results for every root node in vector instead of immediate output
 * for j in range(0, MyFminer.GetNoRootNodes()-1):
 *    result = MyFminer.MineRoot(j);
 *    for i in range(0, result.size()-1):
 *        print result[i];
 * \endcode
 *
 * \subsubsection Ruby Ruby
 *
 * This example assumes that you have created ruby bindings using <code>make ruby</code>.
 * \code
 *
 * require 'last'
 * MyFminer = Last::Last.new()
 * MyFminer.SetMaxHops(25)
 * MyFminer.SetConsoleOut(false)
 * # Add compounds below. IMPORTANT! Do not change settings after adding compounds!
 * MyFminer.AddCompound("O=C(C(C(C=C3)=CC=C3O)=CO2)C1=C2C=C(O)C=C1O" , 1)
 * MyFminer.AddCompound("Oc1ccc(cc1)[C@@H]2Cc3ccc(O)cc3OC2" , 2)
 * MyFminer.AddCompound("O=C1C(C3=CC=C(O)C=C3)=COC2=C1C=CC(O)=C2" , 3)
 * MyFminer.AddCompound("O=C1C(C3=CC=C(OC)C=C3)=COC2=C1C=CC(O)=C2" , 4)
 * MyFminer.AddCompound("OC1=CC=C(CCCCCCCC)C=C1" , 5)
 * MyFminer.AddCompound("C1(C=CC=CC=1C(=C(Cl)Cl)C2=CC=C(C=C2)Cl)Cl" , 6)
 * MyFminer.AddCompound("O=C(C1=C(C=CC=C1)C(=O)OCC(CCCC)CC)OCC(CCCC)CC" , 7)
 * MyFminer.AddCompound("Oc1cc(O)cc2CCCCC[C@@H](O)CCC[C@H](C)OC(=O)c12" , 8)
 * MyFminer.AddCompound("O=C1C2=C(C=C(C=C2O)O)OC(=C1O)C3=CC(=C(C=C3)O)O" , 9)
 * MyFminer.AddCompound("C1(=C(C(=O)C2=C(O1)C=C(C=C2)O)O)C3=CC(O)=C(C=C3)O" , 10)
 * # ... continue adding compounds
 * MyFminer.AddActivity(1.0, 1)
 * MyFminer.AddActivity(1.0, 2)
 * MyFminer.AddActivity(1.0, 3)
 * MyFminer.AddActivity(1.0, 4)
 * MyFminer.AddActivity(1.0, 5)
 * MyFminer.AddActivity(1.0, 6)
 * MyFminer.AddActivity(1.0, 7)
 * MyFminer.AddActivity(1.0, 8)
 * MyFminer.AddActivity(0.0, 9)
 * MyFminer.AddActivity(0.0, 10)
 * # ... continue adding activities (true for active, false for inactive)
 * print MyFminer.GetNoCompounds()  
 * puts " compounds"
 * # gather results for every root node in vector instead of immediate output
 * (0 .. MyFminer.GetNoRootNodes()-1).each do |j|
 *    result = MyFminer.MineRoot(j)
 *    puts "Results"
 *    result.each do |res|
 *        puts res
 *   end
 * end
 *
 * \endcode
 *
 * \subsubsection Const Description of Constructors and Options
 * 
 * For the purpose of demonstration we used a toy database of two compounds and an unusual parameter configuration. Please note, that in general the defaults set by the standard constructor are sensible for most databases. They switch on LAST-PM for 95% significance and a minimum frequency of 2. The complete standard settings are:
 *
 * <ul>
 *  <li>Minimum frequency: <b>2</b></li>
 *  <li>Feature type: <b>Trees</b></li>
 *  <li>Console output: <b>true</b></li>
 *  <li>Aromatic perception: <b>true</b></li>
 *  <li>Refine Singles: <b>false</b></li>
 *  <li>Do Output: <b>true</b></li>
 *  <li>Maximum hops: <b>25</b></li>
 * </ul>
 *
 * \code
 *  //! Constructor for standard settings: 95% significance Lastlevel, minimum frequency 2, type trees
 *  Last ();
 *  \endcode

 * It is recommended to set the maximum number of hops as a first step when too many features are generated or calculation takes too long.
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
 **/
