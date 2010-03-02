/** @mainpage LibLast (liblast)
 *
 * LibLast
 *
 * This is the Last library, available at http://github.com/amaunz/fminer2/tree/master, subdirectory <code>liblast</code>.<br>
 * The Fminer frontend application is available from http://github.com/amaunz/fminer2/tree/master, subdirectory <code>fminer</code>.<br>
 * You may download the scientific documentation from http://cs.maunz.de . The paper is entitled "Large Scale Graph Mining using Backbone Refinement Classes".<br>
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
 *  - OpenBabel: follow the installation instrucations at http://openbabel.org/wiki/Install_(source_code) to build yourself after doing:
 *    \code
 *    apt-get build-dep libopenbabel-dev          # build dependencies for OB
 *    apt-get source libopenbabel-dev             # extracts OB source to the current dir
 *    \endcode
 *    or try the repository version:
 *    \code
 *    apt-get install libopenbabel-dev            # OB binary lib and headers
 *    \endcode
 *  - <a href="http://github.com/amaunz/libfminer/tree" target="_blank">Download the library source code</a> (with git: <code>git clone git://github.com/amaunz/fminer2.git</code>) and cd to <code>liblast</code> subdirectory. Use <code>./configure</code> to configure the Makefile automatically or, in the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - <a href="http://github.com/amaunz/fminer/tree" target="_blank">Download the frontend code</a> (with git: <code>git clone git://github.com/amaunz/fminer2.git</code>) and cd to <code>fminer</code> subdirectory.. Use <code>./configure</code> to configure the Makefile automatically or, in the <code>Makefile</code>, adjust the include (-I) and linker (-L) flags. Run <code>make</code>.
 *  - To create this documentation with doxygen, type 'make doc'. The documentation explains API, constructor usage and options.
 *  @subsection ssec23 Language Portability
 *  The API can be made available to other languages. Follow the installation instructions above. A config file for Swig to automagically create languages bindings exists (<code>rlast_wrap.i</code>). 
 *
 *  The Makefile features a target that creates <b>ruby</b> bindings using this file. On Ubuntu, you can e.g. do this:
 *  - Use <code>./configure <version></code> to configure the Makefile automatically or, adjust the include flags (-I) in the Makefile in the line <code>INCLUDE_RB = ...</code> so that the directory contains file <code>ruby.h</code>. Also, let <code>RUBY = ...</code> point to the right ruby executable.
 *  - Run <code>make ruby</code>.
 *
 *  <a name="Guidance">
 * @section Guidance Guidance on Using (Lib)Last
 *
 * Most setting are sensible by default, see description of constructors and objects below. 
 *
 * I would suggest to manipulate the minimum frequency only at first. The number of fragments output should not be more than a multitude of the number of input graphs.
 * For most chemical databases, a minimum frequency threshold of 4%-6% will deliver good results. LibLast does not support percentage values, you will have to calculate absolute numbers.
 *
 *  @subsection sec3 Examples using the LibLast API
 *  LibLast uses the 'singleton' design pattern known from software engineering, i.e., class instantiation is restricted to one object. To empty the database after a run to feed new compounds, use the Last::Reset() routine. 
 *
 *  The following code demonstrate the use of the Last API from C++ and ruby. It feeds a set of class-labelled molecules in SMILES format (the API currently allows no gSpan input, use the frontend application for that) and calculates a vector of fragments along with statistical relevance and occurrences and prints them out. Every root node corresponds to a single chemical element. The output consists of <a href="http://graphml.graphdrawing.org">GraphML</a> which can be postprocessed to SMARTS patterns using the <a href="http://github.com/amaunz/last-utils" target="_blank">LAST-UTILS</a>.
 *
 * \subsubsection CPP C++
 *
 * This example uses libLast in a C++ program.
 * The example assumes that you have created the C++ library using <code>make</code>.
 *
 * \code
 * #include "last.h"
 * #include <iostream>
 * #include <string.h>
 * using namespace std;
 *
 * Last* MyFminer;
 * int main(int argc, char *argv[], char *envp) {
 *   MyFminer= new Last();
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
 * \subsubsection Ruby Ruby
 *
 * This example assumes that you have created ruby bindings using <code>make ruby</code>.
 * \code
 *
 * require 'last'
 * MyFminer = Last::Last.new()
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
 *
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
 *  <li>Aromatic perception: <b>false</b></li>
 *  <li>Refine Singles: <b>false</b></li>
 *  <li>Do Output: <b>true</b></li>
 *  <li>Separate BBRCs in output by blank line/vector: <b>false</b></li>
 * </ul>
 *
 * \code
 *  //! Constructor for standard settings: 95% significance level, minimum frequency 2, type trees, dynamic upper bound, BBRC
 *  Last ();
 *  \endcode

 * There also exist more flexible constructors:
 * \code
 * //! Like standard constructor, but type and minimum frequency configurable
 * Last (int type, unsigned int minfreq);
 * //! Like standard constructor, but type, minimum frequency, significance level and BBRC configurable
 * Last (int type, unsigned int minfreq, float chisq_val, bool do_backbone);
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
 **/
