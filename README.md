<h1>FAIT</h1>
<b>F</b>FAS based <b>A</b>periodicity detection using <b>I</b>mage processing <b>T</b>echniques

<h2>Dependencies</h2>
You will have to install the following packages:<br/>
<ul>
	<li>Python <a href="www.python.org" target="blank">link</a></li>
	<li>Biopython <a href="http://biopython.org" target="blank">link</a></li>
	<li>Numpy <a href="http://www.numpy.org" target="blank">link</a></li>
	<li>Scipy <a href="http://www.scipy.org" target="blank">link</a></li>
	<li>Matplotlib <a href="http://matplotlib.org/" target="blank">link</a></li>
</ul>
These dependencies are typically available in your standard package management system such as <i>apt</i> or <i>yum</i> on Linux or <i>macports</i> on OSX.
Windows is currently not supported.

<h2>Installation</h2>
Clone / download this repositoy. You will also have to download FFAS from <a href="http://ffas.sanfordburnham.org/ffas-cgi/cgi/download.pl?ses=&rv=&lv=" taget="blank">here</a> and follow the respective installation procedure.

<h2>Running FAIT</h2>
<ol>
	<li>Run FFAS to generate the profile-profile scoring matrix by running <code>ffasMatrix.sh reference.fasta query.fasta pathToOutputMatrix pathToOutputAlignment</code></li>
	<li>Run FAIT to process the profile-profile scoring matrix by <code>fait.py -f pathToOutputMatrix -q query.fasta </code> </li>
</ol>


<h2>Usage</h2>
