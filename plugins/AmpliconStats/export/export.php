<?php
if( isset($_REQUEST["exportdata"]) ) {
	header("Cache-Control: must-revalidate, no-store, post-check=0, pre-check=0");
	header("Content-Length: ".strlen($_REQUEST["exportdata"]));
	header("content-type: plain/text");
	header("content-Disposition: attachment; filename=".$_REQUEST["exportfn"]);
	header("Pragma: no-cache");
	echo $_REQUEST["exportdata"];
}
?>
