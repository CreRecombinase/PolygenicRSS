#!/gpfs/data/xhe-lab/nwk/miniconda3/envs/rssp/bin/perl

use PerlIO::gzip;
use warnings;


open FOO, "<:gzip", $ARGV[0] or die $!;
$outputPrefix = $ARGV[1];
print "output prefix is ${outputPrefix}\n";



$lineNo= $ARGV[2];
print "split every  $lineNo lines\n";

$fileno = 1;
$line = 0;



while (<FOO>) {
    if (!$fh || $line >= $lineNo) {
        open $fh, '>:gzip', "${outputPrefix}_$fileno.gz";
        $fileno++;
        $line = 0;
    }
    print $fh $_; $line++;
}
