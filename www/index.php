
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>	    
	<link href="book/css/bootstrap.min.css" rel="stylesheet">	
        <link href="book/css/R.css" rel="stylesheet">			

  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>

<!-- own website starts here -->

     <br />	 	 
     <h1 align="center">The R-Forge R package copula project</h1> 
     <br />   
    
    <p>
Details and documentation about the <a href="https://cran.r-project.org/package=copula">R package copula</a> can be found on its webpage. The main features of the package are also described in the book <a href="http://copula.r-forge.r-project.org/book/">Elements of Copula Modeling with R</a>.
    </p>      

</body>	      
</html>
