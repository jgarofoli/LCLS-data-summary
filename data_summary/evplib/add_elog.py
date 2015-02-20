import logging
import event_process

class add_elog(event_process.event_process):
    def __init__(self):
        self.output = event_process.event_process_output()
        self.reducer_rank = 0
        self.logger = logging.getLogger(__name__+'.add_available_data')
        return
    
    def beginJob(self):
        if self.parent.rank == self.reducer_rank:
            self.expNum = self.parent.ds.env().expNum()
        return

    def endJob(self):
        if self.parent.rank == self.reducer_rank:
            self.output['in_report']= 'meta'
            self.output['in_report_title'] = 'Elog'
            self.output['text'] = ["<a href=https://pswww.slac.stanford.edu/apps/portal/index.php?exper_id={:0.0f}>Elog</a>".format(self.expNum),]
            self.parent.output.append(self.output)
        return

