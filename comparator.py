"""
Classes for use when comparing functions/graphs

Example usage:

>>> A = Contribution("a.root", "li1corr_eta_0_0.348")
>>> B = Contribution("b.root", "li1corr_eta_0_0.348")
>>> p = Plot([A, B], "graphfunction", xlim=[0, 50])
>>> p.plot()
>>> p.save("AB.pdf")

"""


import os
import numpy as np
import ROOT
import common_utils as cu
from array import array
from MyStyle import My_Style
from collections import OrderedDict


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
# ROOT.gStyle.SetPalette(55)
My_Style.cd()


class MultiFunc(object):
    """Class to handle multiple TF1s, like a TMultiGraph, because ROOT doesn't.

    NOT SOPHISTICATED, JUST A CONTAINER.
    """
    def __init__(self):
        self.funcs = []

    def Add(self, func):
        self.funcs.append(func)

    def Draw(self, option=None):
        """Draw all the TF1s, adjust the axis sizes properly"""
        # if ROOT.gPad:
        #     if not ROOT.gPad.IsEditable():
        #         ROOT.gROOT.MakeDefCanvas()
        #     if "a" in option.lower():
        #         ROOT.gPad.Clear()
        # Draw, and figure out max/min of axes for some auto-ranging
        x_low, x_high = ROOT.Double(0), ROOT.Double(0)
        x_min, x_max = 999, -999
        y_min, y_max = 999, -999
        option = option or ""
        for i, f in enumerate(self.funcs):
            if i == 0:
                f.Draw(option)
            else:
                f.Draw(option + "SAME")

            f.GetRange(x_low, x_high)
            x_min = x_low if x_low < x_min else x_min
            x_max = x_high if x_high > x_max else x_max
            y_min = f.GetMinimum(x_low, x_high) if f.GetMinimum(x_low, x_high) < y_min else y_min
            y_max = f.GetMaximum(x_low, x_high) if f.GetMaximum(x_low, x_high) > y_max else y_max
        if x_max < x_min:
            raise Exception("MultiFunc: x_min > x_max")
        if y_max < y_min:
            raise Exception("MultiFunc: y_min > y_max")

        print (x_min, x_max)
        print (y_min, y_max)

        self.funcs[0].GetXaxis().SetRangeUser(x_min * 0.95, x_max * 1.05)
        self.funcs[0].GetYaxis().SetRangeUser(y_min * 0.95, y_max * 1.05)

    def Mod(self):
        """If we want to do any sort of styling, need to access the first drawn func"""
        return self.funcs[0]


def grab_obj(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    # TODO: checks!
    input_file = cu.open_root_file(file_name)
    obj = cu.get_from_file(input_file, obj_name)
    if isinstance(obj, (ROOT.TH1, ROOT.TGraph)):
        obj.SetDirectory(0)  # Ownership kludge
        input_file.Close()
        return obj.Clone(ROOT.TUUID().AsString())
    else:
        return obj

class Contribution(object):
    """Basic class to handle information about one contribution to a canvas."""

    def __init__(self, obj, label=None,
                 line_width=1, line_color=ROOT.kRed, line_style=1,
                 fill_color=ROOT.kRed, fill_style=0,
                 marker_size=1, marker_color=ROOT.kRed, marker_style=1,
                 normalise_hist=False, rebin_hist=None):
        """
        obj : TH1, TGraph, ...
            Object to be plotted
        label: str
            Title of contribution, to be used in legend.
        line_width: int
            Width of line.
        line_color: int or ROOT color
            Color of line.
        line_style: int
            Line style.
        fill_color: int or ROOT color
            Color of fill.
        fill_style: int
            Fill style.
        marker_size: int
            Size of marker
        marker_color: int or ROOT color
            Color of markers.
        marker_style: int
            Marker style.
        normalise_hist : bool
            If a histogram, normalise so integral = 1
        rebin_hist : int, None
            If a histogram, specify the number of bins to be grouped together.
        """
        self.obj = obj
        self.label = label or ""
        self.line_width = line_width
        self.line_color = line_color
        self.line_style = line_style
        self.fill_color = fill_color
        self.fill_style = fill_style
        self.marker_size = marker_size
        self.marker_color = marker_color
        self.marker_style = marker_style

        self.obj.SetLineWidth(self.line_width)
        self.obj.SetLineColor(self.line_color)
        self.obj.SetLineStyle(self.line_style)
        self.obj.SetFillColor(self.fill_color)
        self.obj.SetFillStyle(self.fill_style)
        self.obj.SetMarkerSize(self.marker_size)
        self.obj.SetMarkerColor(self.marker_color)
        self.obj.SetMarkerStyle(self.marker_style)

        # Match fit to hist styles
        if self.obj.GetListOfFunctions().GetSize() == 1:
            func = self.obj.GetListOfFunctions().At(0)
            func.SetLineStyle(line_style)
            func.SetLineWidth(line_width)
            func.SetLineColor(line_color)

        if rebin_hist and rebin_hist != 1:
            self.obj.Rebin(rebin_hist) # Does this handle 2D hists?
        if normalise_hist and obj.Integral() != 0:
            self.obj.Scale(1./(obj.GetBinWidth(1) * obj.Integral()))
        if isinstance(self.obj, ROOT.TH1) or isinstance(self.obj, ROOT.TH2):
            self.obj.SetDirectory(0)

    def __eq__(self, other):
        return self.obj == other.obj


class Plot(object):
    """
    Basic class to handle information about one plot,
    which can have several contributions.
    """

    def __init__(self, contributions=None, what="graph",
                 title=None, xtitle=None, ytitle=None, xlim=None, ylim=None,
                 legend=True, extend=False,
                 subplot=None, subplot_type=None, subplot_title=None):
        """
        contributions: list
            List of Contribution objects.
        what: str
            Options are "graph", "function", "both".
        title: str
            Title of plot.
        xtitle: str
            X axis title.
        ytitle: str
            Y axis title.
        xlim: list
            Limits of x axis. If None then determines suitable limits.
        ylim: list
            Limits of y axis. If None then determines suitable limits.
        legend: bool
            Include legend on plot.
        extend: bool
            Extend functions to cover the whole x axis range.
        subplot : object, None
            If not None, draws a plot below the main hist with
            all plots compared to the object provided as the argument here.
        subplot_type : str
            The method of comparison in the subplot: ratio, difference, or ddelta/dlambda
        """
        self.contributions = contributions if contributions else []
        self.contributions_objs = []
        options = ['hist', 'graph', 'function', 'both']
        if what not in options:
            raise RuntimeError("`what` argument must be one of %s" % options)
        self.plot_what = what
        self.title = title or ""
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.xlim = xlim
        self.ylim = ylim
        self.do_legend = legend
        self.legend = ROOT.TLegend(0.65, 0.6, 0.94, 0.85) if legend else None
        self.do_extend = extend
        self.container = None
        self.canvas = None
        self.default_canvas_size = (800, 800)
        self.main_pad = None
        self.subplot = subplot
        if subplot_type and subplot_type not in ['ratio', 'diff', "ddelta"]:
            raise RuntimeError("subplot_type must be one of None, ratio, diff, or ddelta")
        self.subplot_type = subplot_type
        if not self.subplot_type:
            self.subplot = None
        self.subplot_container = None
        self.subplot_contributions = []
        self.subplot_pad = None
        self.subplot_ratio_lim = (0, 2)
        self.subplot_pad_height = 0.32
        self.subplot_pad_fudge = 0.01  # to get non-overlapping subplot axis
        self.subplot_line = None  # need this to remain visible...
        self.subplot_title = subplot_title

    def add_contribution(self, *contribution):
        """Add Contribution to Plot. Can be single item or list."""
        self.contributions.extend(*contribution)

    def _create_container(self):
        """Create a container object for lots of the same TObject"""
        if self.plot_what in ["graph", "both"]:
            self.container = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
            self.subplot_container = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
        elif self.plot_what == "function":
            self.container = MultiFunc()
            self.subplot_container = MultiFunc()
        elif self.plot_what == "hist":
            self.container = ROOT.THStack(ROOT.TUUID().AsString(), "")
            self.subplot_container = ROOT.THStack(ROOT.TUUID().AsString(), "")

        if self.container:
            ROOT.SetOwnership(self.container, False)
        if self.subplot_container:
            ROOT.SetOwnership(self.subplot_container, False)


    def _populate_container_and_legend(self):
        """Add objects to the container, and to the legend"""
        ref_func = None
        if (self.subplot and (self.plot_what == "both" or self.plot_what == "function")):
            ref_func = self.subplot.obj.GetListOfFunctions().At(0)

        ref_graph = None
        if (self.subplot and (self.plot_what == "both" or self.plot_what == "graph")):
            ref_graph = self.subplot.obj

        for contrib in self.contributions:
            if self.plot_what == "graph":
                # if drawing only graph, we need to remove the fit if there is one
                if (contrib.obj.GetListOfFunctions().GetSize() > 0):
                    contrib.obj.GetListOfFunctions().Remove(contrib.obj.GetListOfFunctions().At(0))
            else:
                # if drawing function, extend range as per user's request
                if self.do_extend:
                    if self.plot_what == 'function':
                        contrib.obj.SetRange(self.xlim[0], self.xlim[1])
                    elif self.plot_what == "both":
                        contrib.obj.GetListOfFunctions().At(0).SetRange(self.xlim[0], self.xlim[1])
            self.container.Add(contrib.obj)
            self.contributions_objs.append(contrib.obj)

            if self.do_legend:
                # Split text by newline \n
                # Add an entry for each line
                opt = "LP" if self.plot_what == "hist" else "LP"
                for i, substr in enumerate(contrib.label.split("\n")):
                    if i == 0:
                        self.legend.AddEntry(contrib.obj, substr, opt)
                    else:
                        self.legend.AddEntry(0, substr, "")


            # Add contributions for the subplot
            if self.subplot:
                subplot_obj = self.subplot.obj.Clone()
                if contrib != self.subplot:
                    new_obj = contrib.obj.Clone()
                    if (self.subplot_type == "ratio"):
                        if self.plot_what == 'graph' or self.plot_what == 'both':
                            interp_obj = self._construct_interpolated_graph(contrib.obj, ref_graph)
                            new_obj = self._construct_divided_graph(interp_obj, ref_graph)
                            new_obj.SetLineWidth(contrib.line_width)
                            new_obj.SetLineColor(contrib.line_color)
                            new_obj.SetLineStyle(contrib.line_style)
                            new_obj.SetFillColor(contrib.fill_color)
                            new_obj.SetFillStyle(contrib.fill_style)
                            new_obj.SetMarkerSize(contrib.marker_size)
                            new_obj.SetMarkerColor(contrib.marker_color)
                            new_obj.SetMarkerStyle(contrib.marker_style)
                            if self.plot_what == "both":
                                orig_func = contrib.obj.GetListOfFunctions().At(0)
                                new_func = self._construct_divided_function(orig_func, ref_func)
                                new_func.SetLineColor(orig_func.GetLineColor())
                                new_func.SetLineWidth(orig_func.GetLineWidth())
                                new_func.SetLineStyle(orig_func.GetLineStyle())

                        elif self.plot_what == 'hist':
                            new_obj.Divide(subplot_obj)

                    elif (self.subplot_type == "diff"):
                        if self.plot_what == 'graph':
                            pass
                        elif self.plot_what == 'hist':
                            new_obj.Add(subplot_obj, -1.)

                    elif (self.subplot_type == "ddelta"):
                        # Do the differntial delta spectrum, see 1704.03878
                        new_obj.Add(subplot_obj, -1.)
                        new_obj.Multiply(new_obj)
                        sum_hist = contrib.obj.Clone()
                        sum_hist.Add(subplot_obj)
                        new_obj.Divide(sum_hist)
                        new_obj.Scale(0.5)

                    self.subplot_container.Add(new_obj)
                    self.subplot_contributions.append(new_obj)

    def _construct_interpolated_graph(self, from_this_graph, ref_graph):
        """Construct a new graph from from_this_graph using the x values from ref_graph as reference + interpolation"""
        ref_x, ref_y = cu.get_xy(ref_graph)
        this_x, this_y = cu.get_xy(from_this_graph)
        this_ex, this_ey = cu.get_exey(from_this_graph)
        new_x, new_y = [], []
        new_ex, new_ey = [], []
        for x in ref_x:

            # ignore anything outside the original range of from_this_graph
            if x > max(this_x) or x < min(this_x):
                continue

            tmp_y = from_this_graph.Eval(x)

            new_x.append(x)
            new_y.append(tmp_y)
            new_ex.append(0.)
            new_ey.append(0.)

        new_graph = ROOT.TGraphErrors(len(new_x), array('d', new_x), array('d', new_y))
        return new_graph

    def _construct_divided_graph(self, this_graph, ref_graph):
        """Construct graph by doing this_graph/ref_graph"""
        one_dict = {k:v for k, v in zip(*cu.get_xy(this_graph))}
        one_dict = OrderedDict(sorted(one_dict.items(), key=lambda t: t[0]))

        other_dict = {k:v for k, v in zip(*cu.get_xy(ref_graph))}
        other_dict = OrderedDict(sorted(other_dict.items(), key=lambda t: t[0]))

        new_x, new_y = [], []
        for x, y in one_dict.items():
            if x in other_dict:
                new_x.append(x)
                new_y.append(y / other_dict[x])

        gr = ROOT.TGraphErrors(len(new_x), array('f', new_x), array('f', new_y))
        return gr

    def _construct_divided_function(self, this_func, ref_func):
        """Construct function by doing this_func / ref_func"""
        this_formula = this_func.GetExpFormula("CLINGP")
        ref_formula = this_func.GetExpFormula("CLINGP")
        new_formula = "(%s)/(%s)" % (this_formula, ref_formula)
        new_func = TFormula("new_formula" + ROOT.TUUID().AsString(), new_formula)
        return new_func

    def _style_legend(self):
        # self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)

    def _rescale_plot_labels(self, container, factor):
        # What a pile of wank, why does ROOT scale all these sizes?
        container.GetXaxis().SetLabelSize(container.GetXaxis().GetLabelSize()/factor)
        container.GetXaxis().SetTitleSize(container.GetXaxis().GetTitleSize()/factor)
        container.GetXaxis().SetTitleOffset(container.GetXaxis().GetTitleOffset()*factor)  # doesn't seem to work?
        container.GetXaxis().SetTickLength(container.GetXaxis().GetTickLength()/factor)

        container.GetYaxis().SetLabelSize(container.GetYaxis().GetLabelSize()/factor)
        container.GetYaxis().SetTitleSize(container.GetYaxis().GetTitleSize()/factor)
        container.GetYaxis().SetTitleOffset(container.GetYaxis().GetTitleOffset()*factor)
        # container.GetYaxis().SetTickLength(0.03/factor)

    def set_logx(self, state=True):
        self.main_pad.SetLogx(int(state))
        self.subplot_pad.SetLogx(int(state))
        if self.container:
            ax = self.container.GetXaxis()
            if ax:
                ax.SetMoreLogLabels()
        if self.subplot_container:
            ax = self.subplot_container.GetXaxis()
            if ax:
                ax.SetMoreLogLabels()

    def set_logy(self, state=True):
        self.main_pad.SetLogy(int(state))
        self.subplot_pad.SetLogy(int(state))
        if self.container:
            ax = self.container.GetYaxis()
            if ax:
                ax.SetMoreLogLabels()
        if self.subplot_container:
            ax = self.subplot_container.GetYaxis()
            if ax:
                ax.SetMoreLogLabels()

    def set_logz(self, state=True):
        self.main_pad.SetLogz(int(state))

    def plot(self, draw_opts=None):
        """Make the plot.

        draw_opts: str
            Same options as you would pass to Draw() in ROOT.
        """
        if not self.container:
            # First make a container
            self._create_container()

            # Now add all the contributions to the container, styling as we go
            if len(self.contributions) == 0:
                raise UnboundLocalError("contributions list is empty")

            self._populate_container_and_legend()

        # Set default drawing opts
        # need different default drawing options for TF1s vs TGraphs
        # (ALP won't actually display axis for TF1)
        if draw_opts is None:
            if self.plot_what in ["graph", 'both']:
                draw_opts = "ALP"
            elif self.plot_what in ['function']:
                draw_opts = ""
            elif self.plot_what == "hist":
                draw_opts = ""

        # Need a canvas before drawing
        # If we have "SAME" then we want to draw this Plot object
        # on the current canvas
        # Otherwise, we create a new one
        if not self.canvas:
            if "SAME" in draw_opts:
                self.canvas = ROOT.gPad
                self.canvas.cd()
                print ("Using existing canvas", self.canvas.GetName())
            else:
                self.canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", *self.default_canvas_size)
                self.canvas.SetTicks(1, 1)
                right_margin = 0.03
                if self.subplot:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, self.subplot_pad_height+self.subplot_pad_fudge, 1, 1)
                    ROOT.SetOwnership(self.main_pad, False)
                    self.main_pad.SetTicks(1, 1)
                    self.main_pad.SetBottomMargin(2*self.subplot_pad_fudge)
                    self.main_pad.SetTopMargin(0.12)
                    self.main_pad.SetRightMargin(right_margin)
                    self.canvas.cd()
                    self.main_pad.Draw()
                    self.subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, self.subplot_pad_height-self.subplot_pad_fudge)
                    ROOT.SetOwnership(self.subplot_pad, False)
                    self.subplot_pad.SetTicks(1, 1)
                    self.subplot_pad.SetFillColor(0)
                    self.subplot_pad.SetFillStyle(0)
                    self.subplot_pad.SetTopMargin(2*self.subplot_pad_fudge)
                    self.subplot_pad.SetRightMargin(right_margin)
                    self.subplot_pad.SetBottomMargin(0.35)
                    self.canvas.cd()
                    self.subplot_pad.Draw()
                else:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, 0, 1, 1)
                    ROOT.SetOwnership(self.main_pad, False)
                    self.main_pad.SetRightMargin(right_margin)
                    self.main_pad.SetTicks(1, 1)
                    self.main_pad.Draw()

        self.canvas.cd()
        self.main_pad.cd()

        # Plot the container
        self.container.Draw(draw_opts)

        # Customise
        if self.plot_what != 'function':
            modifier = self.container
        else:
            modifier = self.container.Mod()

        # Use the x/y axis labels from the first contribution as defaults
        if not self.xtitle:
            self.xtitle = self.contributions_objs[0].GetXaxis().GetTitle()
        if not self.ytitle:
            self.ytitle = self.contributions_objs[0].GetYaxis().GetTitle()

        modifier.SetTitle("%s;%s;%s" % (self.title, self.xtitle, self.ytitle))

        for ob in self.contributions_objs:
            ob.GetXaxis().SetTitle(self.xtitle)
            ob.GetXaxis().SetTitle(self.ytitle)

        if self.xlim:
            if self.plot_what == "graph":
                modifier.GetXaxis().SetLimits(*self.xlim)
            else:
                modifier.GetXaxis().SetRangeUser(*self.xlim)
        if self.ylim:
            # dont use the SetLimits for graphs, that doesnt work properly.
            # no idea, ROOT is fucking insane
            modifier.GetYaxis().SetRangeUser(*self.ylim)
            modifier.SetMinimum(self.ylim[0])  # use this, not SetRangeUser()
            modifier.SetMaximum(self.ylim[1])

        if self.main_pad.GetLogx():
            modifier.GetXaxis().SetMoreLogLabels()
        if self.main_pad.GetLogy():
            modifier.GetYaxis().SetMoreLogLabels()

        # Draw it again to update
        if self.plot_what == "graph":
            self.container.Draw(draw_opts)

        # Plot legend
        if self.do_legend:
            self._style_legend()
            self.canvas.cd()
            self.legend.Draw()

        if self.subplot:
            self._rescale_plot_labels(modifier, 1-self.subplot_pad_height)

            # Get rid of main plot x axis labels
            modifier.GetHistogram().GetXaxis().SetLabelSize(0)
            modifier.GetXaxis().SetLabelSize(0)
            modifier.GetHistogram().GetXaxis().SetTitleSize(0)
            modifier.GetXaxis().SetTitleSize(0)

            self.subplot_pad.cd()
            self.subplot_container.Draw(draw_opts)

            if self.subplot_title == None:
                if (self.subplot_type == "ratio"):
                    self.subplot_title = "#splitline{Ratio vs}{%s}" % (self.subplot.label)
                elif (self.subplot_type == "diff"):
                    self.subplot_title = "#splitline{Difference}{vs %s}" % (self.subplot.label)
                elif (self.subplot_type == "ddelta"):
                    self.subplot_title = "d#Delta/d#lambda"

            self.subplot_container.SetTitle(";%s;%s" % (self.xtitle, self.subplot_title))

            self._rescale_plot_labels(self.subplot_container, self.subplot_pad_height)
            self.subplot_container.GetXaxis().SetTitleOffset(self.subplot_container.GetXaxis().GetTitleOffset()*2.8)
            self.subplot_container.GetYaxis().SetNdivisions(505)

            # sync x limits with main plot
            main_xlimits = [self.container.GetXaxis().GetXmin(), self.container.GetXaxis().GetXmax()]
            self.subplot_container.GetXaxis().SetLimits(*main_xlimits)  # setrangeuser doesn't work for x axis

            if self.subplot_pad.GetLogx():
                self.subplot_container.GetXaxis().SetMoreLogLabels()
            if self.subplot_pad.GetLogy():
                self.subplot_container.GetYaxis().SetMoreLogLabels()

            if self.subplot_type == "ratio":

                if self.plot_what == "hist":
                    # self.subplot_container.SetMinimum(self.subplot_ratio_lim[0])  # use this, not SetRangeUser()
                    self.subplot_container.SetMinimum(0)  # use this, not SetRangeUser()

                    # Make sure that the upper limit is the largest bin of the contributions,
                    # so long as it is within 1.5 and some upper limit
                    bin_meds = [np.max(cu.th1_to_arr(h)) for h in self.subplot_contributions]
                    self.subplot_container.SetMaximum(min(10, max(1.5, 1.*max(bin_meds))))

                    # Make sure the lower limit is the smallest bin of the contributions,
                    # so long as it is within 0 and 0.5
                    bin_mins = [np.min(cu.th1_to_arr(h)) for h in self.subplot_contributions]
                    self.subplot_container.SetMinimum(min(0.5, min(bin_mins)))

                    # self.subplot_container.SetMaximum(1.5)
                    # self.subplot_container.SetMinimum(0.5)

                xax = modifier.GetXaxis()
                self.subplot_line = ROOT.TLine(xax.GetXmin(), 1., xax.GetXmax(), 1.)
                self.subplot_line.SetLineStyle(2)
                self.subplot_line.SetLineWidth(2)
                self.subplot_line.SetLineColor(ROOT.kBlack)
                self.subplot_line.Draw()

                self.subplot_line_up = ROOT.TLine(xax.GetXmin(), 1.02, xax.GetXmax(), 1.02)
                self.subplot_line_up.SetLineStyle(3)
                self.subplot_line_up.SetLineWidth(2)
                self.subplot_line_up.SetLineColor(ROOT.kBlack)
                self.subplot_line_up.Draw()

                self.subplot_line_down = ROOT.TLine(xax.GetXmin(), 0.98, xax.GetXmax(), 0.98)
                self.subplot_line_down.SetLineStyle(3)
                self.subplot_line_down.SetLineWidth(2)
                self.subplot_line_down.SetLineColor(ROOT.kBlack)
                self.subplot_line_down.Draw()

            self.subplot_pad.Update()
            self.canvas.Update()


    def save(self, filename):
        """Save the plot to file. Do some check to make sure dir exists."""
        filename = os.path.abspath(filename)
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        self.canvas.SaveAs(filename)
        del self.canvas
