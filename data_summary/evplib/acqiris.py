import event_process
import psana
import os
import logging
import toolbox
import pylab
from common import strtype


class acqiris(event_process.event_process):
    def __init__(self):
        self.data = {}
        self.output = {}
        self.logger = logging.getLogger(__name__+'.acqiris')
        return

    def set_stuff(self,src,psana_device,in_report=None,in_report_title=None,histmin=400,histmax=450):
        self.src = psana.Source(src)
        self.dev = psana_device
        self.output['in_report'] = in_report
        self.output['in_report_title'] = in_report_title
        self.histmin = histmin
        self.histmax = histmax
        return

    def replicate_info(self):
        args = ( str(self.src), strtype(self.dev))
        kwargs = { 'in_report': self.output['in_report'], 'in_report_title': self.output['in_report_title'], 
                'histmin': self.histmin, 'histmax': self.histmax }
        self.logger.info('args: {:}'.format(repr(args)))
        self.logger.info('kwargs: {:}'.format(repr(kwargs)))
        return ('set_stuff',args,kwargs)

    def event(self,evt):
        self.raw_traces = evt.get(self.dev,self.src)
        if self.raw_traces is None:
            self.logger.error('No acqiris found in event {:}'.format(self.parent.eventN))
            return
        #self.logger.debug( 'acqiris traces = {:}'.format(self.raw_traces.data_shape() ))
        self.trace = list(self.raw_traces.data(0).waveforms()[0]) # or 5? idk
        self.peak = self.trace.index( max(self.trace) )
        for evr in self.parent.shared['evr']:
            #print "rank {:} evr {:} peak {:}".format(self.parent.rank, evr, peak )
            self.data.setdefault( evr, toolbox.myhist(500,0,500)).fill( self.peak )
        return

    def endJob(self):
        # do the stuff on all the other nodes (e.g. send stuff to reducer_rank)
        self.reduced_data = {}
        for evr in self.data:
            self.logger.info('mpi reducing {:}'.format(evr))
            self.reduced_data[evr] = self.data[evr].reduce(self.parent.comm,self.reducer_rank)


        if self.parent.rank == self.reducer_rank:
            self.output['table'] = {}
            self.output['figures'] = {}
            fig = pylab.figure()
            for evr in sorted(self.reduced_data):
                fig.clear()
                newdata = self.reduced_data[evr]
                self.logger.info( "{:} mean: {:0.2f}, std: {:0.2f}, min {:0.2f}, max {:0.2f}".format( evr, newdata.mean(), newdata.std(), newdata.minval, newdata.maxval ) )
                self.output['table'][evr] = {}
                self.output['table'][evr]['Mean Arrival Bin'] = newdata.mean()
                self.output['table'][evr]['RMS'] = newdata.std()
                self.output['table'][evr]['min'] = newdata.minval
                self.output['table'][evr]['max'] = newdata.maxval

                self.output['figures'][evr] = {}
                #hh = pylab.hist(newdata,bins=100,range=(0,100))
                newX, newY = self.reduced_data[evr].mksteps()
                hh = pylab.fill_between( newX, 0, newY[:-1] )
                pylab.title( 'evr '+repr(evr) )
                pylab.xlim(self.histmin,self.histmax)
                pylab.ylim( 0 , max(self.reduced_data[evr].binentries)*1.1 )
                #pylab.savefig( os.path.join( self.parent.output_dir, 'figure_evr_{:}.pdf'.format( evr ) ) )
                #self.output['figures'][evr]['pdf'] = os.path.join( self.parent.output_dir, 'figure_evr_{:}.pdf'.format( evr ) )
                pylab.savefig( os.path.join( self.parent.output_dir, 'figure_evr_{:}.png'.format( evr ) ) )
                self.output['figures'][evr]['png'] = os.path.join( self.parent.output_dir, 'figure_evr_{:}.png'.format( evr ) )
            del fig
            self.parent.output.append(self.output)
        return
