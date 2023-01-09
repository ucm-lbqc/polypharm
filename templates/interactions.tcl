set folder [lindex $argv 0]
set output [lindex $argv 1]
set pdbfiles [glob $folder/*.mae]
set ligand "UNK"

set resids %{resids}
set radius %{radius}
set chains %{chains}

set outfile [open ${output}_interactions.csv w]

puts -nonewline $outfile "Name"
foreach resid $resids {
    foreach chain $chains {
        puts -nonewline $outfile ",$chain:$resid"
    }
}
puts $outfile ""


foreach file $pdbfiles {
    set basename [file rootname [file tail $file]]
    puts -nonewline $outfile $basename
    mol new $file
    foreach resid $resids {
        foreach chain $chains {
            set sel [atomselect top "same residue as (resid $resids and chain $chain and within $radius of resname $ligand)"]
            set current_resids [lsort -unique [$sel get resid]]
            if {[lsearch $current_resids $resid] > -1} {
                # lset mat $i [expr $j*[llength $chains]+$offset] 1
                puts -nonewline $outfile ",1"
            } else {
                puts -nonewline $outfile ",0"
            }
        }
    }
    puts $outfile ""
    mol delete top
}

close $outfile

exit
