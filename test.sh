        #!/bin/bash
        for i in `seq 1 100`;
        do
                echo $i
                ./bin/DDiMAP -b ~/Google\ Drive/DataExchangeUR/DDiMAP_AsynthExpt_W.bam -f ~/Google\ Drive/DataExchangeUR/refSeqEnhanced.fa -r $i | grep SNV
        done    
