set folder [lindex $argv 0]
set output [lindex $argv 1]

set pdbfiles [glob $folder/*.mae]
set ligand "UNK"

set resids %{resids}
set radius %{radius}

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
        lassign $resid ch num
        set sel [atomselect top "same residue as (resid $num and chain $ch and within $radius of resname $ligand)"]
        if {[$sel num] > 0} {
            puts -nonewline $outfile ",1"
        } else {
            puts -nonewline $outfile ",0"
        }
    }
    puts $outfile ""
    mol delete top
}

close $outfile

exit
