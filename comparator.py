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
from MyStyle import My_Style


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

        self.funcs[0].GetXaxis().SetRangeUser(x_min * 0.95, x_max * 1.05)
        self.funcs[0].GetYaxis().SetRangeUser(y_min * 0.95, y_max * 1.05)

    def Mod(self):
        """If we want to do any sort of styling, need to access the first drawn func"""
        return self.funcs[0]


def grab_obj(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    # TODO: checks!
    input_file = cu.open_root_file(file_name)
    obj = cu.get_from_tfile(input_file, obj_name)
    # print("Getting", obj_name, "from", file_name)
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
                 normalise_hist=False, rebin_hist=None, subplot=None, subplot2=None):
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
        subplot : ROOT.TObject
            Object to use for subplot, if there is one
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
        self.subplot = subplot
        self.subplot2 = subplot2

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
            # self.obj.Scale(1./(obj.GetBinWidth(1) * obj.Integral()))
            self.obj.Scale(1./(obj.Integral()))
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
                 subplot=None, subplot_type=None, subplot_title=None, subplot_limits=None,
                 subplot2=None, subplot_type2=None, subplot_title2=None, subplot_limits2=None,
                 has_data=True):
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
            If you want the subplot for each contrib, this must be None.
        subplot_type : str
            The method of comparison in the subplot: ratio, difference, or ddelta/dlambda
        subplot_limits : (float, float), optional
            Set hard limits on subplot y axis range. If None, will choose some
            vaguely sensible ones
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
        self.legend = ROOT.TLegend(0.65, 0.67, 0.94, 0.92) if legend else None
        self.do_extend = extend
        self.container = None
        self.canvas = None
        # self.default_canvas_size = (600, 800)
        self.default_canvas_size = (800, 600)
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
        self.subplot_limits = subplot_limits
        self.subplot_pad_height = 0.32
        self.subplot_pad_fudge = 0.01  # to get non-overlapping subplot axis
        self.subplot_line = None  # need this to remain visible...
        self.subplot_title = subplot_title

        self.subplot2 = subplot
        if subplot_type2 and subplot_type2 not in ['ratio', 'diff', "ddelta"]:
            raise RuntimeError("subplot_type must be one of None, ratio, diff, or ddelta")
        self.subplot_type2 = subplot_type2
        if not self.subplot_type2:
            self.subplot2 = None
        self.subplot_container2 = None
        self.subplot_contributions2 = []
        self.subplot_pad2 = None
        self.subplot_limits2 = subplot_limits2
        self.subplot_pad_height2 = 0.32
        self.subplot_pad_fudge2 = 0.01  # to get non-overlapping subplot axis
        self.subplot_line2 = None  # need this to remain visible...
        self.subplot_title2 = subplot_title2
        if self.subplot_type2 and self.subplot_type:
            self.subplot_pad_height = 0.12
            self.subplot_pad_height2 = 0.20
        
        self.is_preliminary = True
        self.has_data = has_data

    def add_contribution(self, *contribution):
        """Add Contribution to Plot. Can be single item or list."""
        self.contributions.extend(*contribution)

    def _create_container(self):
        """Create a container object for lots of the same TObject"""
        if self.plot_what in ["graph", "both"]:
            self.container = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
            self.subplot_container = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
            self.subplot_container2 = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
        elif self.plot_what == "function":
            self.container = MultiFunc()
            self.subplot_container = MultiFunc()
            self.subplot_container2 = MultiFunc()
        elif self.plot_what == "hist":
            self.container = ROOT.THStack(ROOT.TUUID().AsString(), "")
            self.subplot_container = ROOT.THStack(ROOT.TUUID().AsString(), "")
            self.subplot_container2 = ROOT.THStack(ROOT.TUUID().AsString(), "")

        if self.container:
            ROOT.SetOwnership(self.container, False)
        if self.subplot_container:
            ROOT.SetOwnership(self.subplot_container, False)
        if self.subplot_container2:
            ROOT.SetOwnership(self.subplot_container2, False)


    def _populate_container_and_legend(self):
        """Add objects to the container, and to the legend"""
        for contrib in self.contributions:
            if self.plot_what == "graph":
                # if drawing only graph, we need to remove the fit if there is one
                if (contrib.obj.GetListOfFunctions().GetSize() > 0):
                    contrib.obj.GetListOfFunctions().Remove(contrib.obj.GetListOfFunctions().At(0))
            else:
                if isinstance(contrib.obj, ROOT.TH1) and contrib.obj.GetEntries() == 0:
                    continue
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
            if self.subplot_type:
                if self.subplot:
                    print(self.subplot)
                    # Use one reference object for all entries
                    subplot_obj = self.subplot.obj.Clone()
                    if contrib != self.subplot:
                        new_hist = contrib.obj.Clone()
                        if (self.subplot_type == "ratio"):
                            new_hist.Divide(subplot_obj)
                        elif (self.subplot_type == "diff"):
                            new_hist.Add(subplot_obj, -1.)
                        elif (self.subplot_type == "ddelta"):
                            # Do the differntial delta spectrum, see 1704.03878
                            new_hist.Add(subplot_obj, -1.)
                            new_hist.Multiply(new_hist)
                            sum_hist = contrib.obj.Clone()
                            sum_hist.Add(subplot_obj)
                            new_hist.Divide(sum_hist)
                            new_hist.Scale(0.5)
                        self.subplot_container.Add(new_hist)
                        self.subplot_contributions.append(new_hist)
                elif contrib.subplot is not None:
                    new_hist = contrib.obj.Clone()
                    ref_obj = contrib.subplot
                    if (self.subplot_type == "ratio"):
                        new_hist.Divide(ref_obj)
                    elif (self.subplot_type == "diff"):
                        new_hist.Add(ref_obj, -1.)
                    elif (self.subplot_type == "ddelta"):
                        # Do the differntial delta spectrum, see 1704.03878
                        new_hist.Add(ref_obj, -1.)
                        new_hist.Multiply(new_hist)
                        sum_hist = contrib.obj.Clone()
                        sum_hist.Add(ref_obj)
                        new_hist.Divide(sum_hist)
                        new_hist.Scale(0.5)
                    self.subplot_container.Add(new_hist)
                    self.subplot_contributions.append(new_hist)
            
            # Add contributions for the subplot2
            if self.subplot_type2:
                if self.subplot2:
                    print(self.subplot2)
                    # Use one reference object for all entries
                    subplot_obj = self.subplot2.obj.Clone()
                    if contrib != self.subplot2:
                        new_hist2 = contrib.obj.Clone()
                        if (self.subplot_type2 == "ratio"):
                            new_hist2.Divide(subplot_obj)
                        elif (self.subplot_type2 == "diff"):
                            new_hist2.Add(subplot_obj, -1.)
                        elif (self.subplot_type2 == "ddelta"):
                            # Do the differntial delta spectrum, see 1704.03878
                            new_hist2.Add(subplot_obj, -1.)
                            new_hist2.Multiply(new_hist2)
                            sum_hist2 = contrib.obj.Clone()
                            sum_hist2.Add(subplot_obj)
                            new_hist2.Divide(sum_hist2)
                            new_hist2.Scale(0.5)
                        self.subplot_container2.Add(new_hist2)
                        self.subplot_contributions2.append(new_hist2)
                elif contrib.subplot2 is not None:
                    new_hist2 = contrib.obj.Clone()
                    ref_obj = contrib.subplot2
                    if (self.subplot_type2 == "ratio"):
                        new_hist2.Divide(ref_obj)
                    elif (self.subplot_type2 == "diff"):
                        new_hist2.Add(ref_obj, -1.)
                    elif (self.subplot_type2 == "ddelta"):
                        # Do the differntial delta spectrum, see 1704.03878
                        new_hist2.Add(ref_obj, -1.)
                        new_hist2.Multiply(new_hist2)
                        sum_hist2 = contrib.obj.Clone()
                        sum_hist2.Add(ref_obj)
                        new_hist2.Divide(sum_hist2)
                        new_hist2.Scale(0.5)
                    self.subplot_container2.Add(new_hist2)
                    self.subplot_contributions2.append(new_hist2)

    def _style_legend(self):
        # self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)

    def _rescale_plot_labels(self, container, factor):
        # Eurgh, why does ROOT scale all these sizes?
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
        if self.subplot_pad:
            self.subplot_pad.SetLogx(int(state))
        if self.subplot_pad2:
            self.subplot_pad2.SetLogx(int(state))
        if self.container:
            ax = self.container.GetXaxis()
            # if ax:
            #     ax.SetMoreLogLabels()
        if self.subplot_container:
            ax = self.subplot_container.GetXaxis()
        if self.subplot_container2:
            ax = self.subplot_container2.GetXaxis()
            # if ax:
            #     ax.SetMoreLogLabels()

    def set_logy(self, state=True):
        self.main_pad.SetLogy(int(state))
        # if self.subplot_pad:
        #     self.subplot_pad.SetLogy(int(state))
        if self.container:
            ax = self.container.GetYaxis()
            # if ax:
            #     ax.SetMoreLogLabels()
        # if self.subplot_container:
        #     ax = self.subplot_container.GetYaxis()
        #     if ax:
        #         ax.SetMoreLogLabels()

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
                draw_opts = "NOSTACK"

        # Need a canvas before drawing
        # If we have "SAME" then we want to draw this Plot object
        # on the current canvas
        # Otherwise, we create a new one
        if not self.canvas:
            if "SAME" in draw_opts:
                self.canvas = ROOT.gPad
                self.canvas.cd()
                # print("Using existing canvas", self.canvas.GetName())
            else:
                self.canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", *self.default_canvas_size)
                self.canvas.SetTicks(1, 1)
                right_margin = 0.03
                top_margin = 0.05
                if self.subplot_type:
                    if not self.subplot_type2:
                        self.main_pad = ROOT.TPad("main_pad", "", 0, self.subplot_pad_height+self.subplot_pad_fudge, 1, 1)
                        ROOT.SetOwnership(self.main_pad, False)
                        self.main_pad.SetTicks(1, 1)
                        self.main_pad.SetBottomMargin(2*self.subplot_pad_fudge)
                        self.main_pad.SetTopMargin(top_margin / (1-self.subplot_pad_height))
                        # self.main_pad.SetRightMargin(right_margin / (1-self.subplot_pad_height))
                        self.main_pad.SetRightMargin(right_margin)
                        self.canvas.cd()
                        self.main_pad.Draw()

                        self.subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, self.subplot_pad_height-self.subplot_pad_fudge)
                        ROOT.SetOwnership(self.subplot_pad, False)
                        self.subplot_pad.SetTicks(1, 1)
                        self.subplot_pad.SetFillColor(0)
                        self.subplot_pad.SetFillStyle(0)
                        self.subplot_pad.SetTopMargin(4*self.subplot_pad_fudge)
                        self.subplot_pad.SetRightMargin(right_margin)
                        self.subplot_pad.SetBottomMargin(0.35)
                        self.canvas.cd()
                        self.subplot_pad.Draw()
                    else:
                        self.main_pad = ROOT.TPad("main_pad", "", 0, self.subplot_pad_height+self.subplot_pad_fudge+self.subplot_pad_height2+self.subplot_pad_fudge2, 1, 1)
                        ROOT.SetOwnership(self.main_pad, False)
                        self.main_pad.SetTicks(1, 1)
                        self.main_pad.SetBottomMargin(2*self.subplot_pad_fudge)
                        self.main_pad.SetTopMargin(top_margin / (1-self.subplot_pad_height-self.subplot_pad_height2))
                        # self.main_pad.SetRightMargin(right_margin / (1-self.subplot_pad_height))
                        self.main_pad.SetRightMargin(right_margin)
                        self.canvas.cd()
                        self.main_pad.Draw()

                        # I have no idea how to do this properly, all hacked until looks good
                        self.subplot_pad = ROOT.TPad("subplot_pad", "", 0, self.subplot_pad_height2, 1, self.subplot_pad_height2+self.subplot_pad_fudge2+self.subplot_pad_height+self.subplot_pad_fudge)
                        ROOT.SetOwnership(self.subplot_pad, False)
                        self.subplot_pad.SetTicks(1, 1)
                        self.subplot_pad.SetFillColor(0)
                        self.subplot_pad.SetFillStyle(0)
                        # self.subplot_pad.SetTopMargin(4*self.subplot_pad_fudge)
                        self.subplot_pad.SetRightMargin(right_margin)
                        # self.subplot_pad.SetBottomMargin(0.35*0.5)
                        self.canvas.cd()
                        self.subplot_pad.Draw()

                        self.subplot_pad2 = ROOT.TPad("subplot_pad2", "", 0, 0, 1, self.subplot_pad_height2+self.subplot_pad_fudge2)
                        ROOT.SetOwnership(self.subplot_pad2, False)
                        self.subplot_pad2.SetTicks(1, 1)
                        self.subplot_pad2.SetFillColor(0)
                        self.subplot_pad2.SetFillStyle(0)
                        # self.subplot_pad2.SetTopMargin(3*self.subplot_pad_fudge2)
                        self.subplot_pad2.SetRightMargin(right_margin)
                        self.subplot_pad2.SetBottomMargin(0.35*1.2)
                        self.canvas.cd()
                        self.subplot_pad2.Draw()
                else:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, 0, 1, 1)
                    ROOT.SetOwnership(self.main_pad, False)
                    self.main_pad.SetRightMargin(right_margin)
                    self.main_pad.SetTopMargin(top_margin)
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

        modifier.SetTitle(";%s;%s" % (self.xtitle, self.ytitle))

        for ob in self.contributions_objs:
            # who knows which actually sets the thing?
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
        else:
            # add some padding
            self.canvas.Update()
            # urgh why doesnt THStack.GetMaximum() return the actual maximum
            # GetYaxis().GetXmax() doesnt work either
            ymax = max([o.GetMaximum() for o in self.contributions_objs])
            modifier.SetMaximum(ymax * 1.35)  # TODO different for log and lin scales

        # Draw it again to update
        if self.plot_what == "graph":
            self.container.Draw(draw_opts)

        # Plot legend
        if self.do_legend:
            self._style_legend()
            self.canvas.cd()
            self.legend.Draw()
            self.main_pad.cd()


        # Add CMS text
        self.canvas.cd()
        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        cms_latex.SetTextSize(0.035)
        latex_height = 0.97
        # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
        if self.is_preliminary:
            if self.has_data:
                cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
            else:
                # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Simulation}")
                cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary Simulation}")
        # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}")
        cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        if self.has_data:
            cms_latex.DrawLatex(0.97, latex_height, " xxx fb^{-1} (13 TeV)")
        else:
            cms_latex.DrawLatex(0.97, latex_height, " (13 TeV)")

        # lumi_latex = ROOT.TLatex()
        # lumi_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        # lumi_latex.SetTextFont(42)
        # lumi_latex.SetTextSize(0.035)
        # lumi_latex.DrawLatex(0.97, 0.93, "35.9 fb^{-1} (13 TeV)")
        # lumi_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        # lumi_latex.DrawLatex(0.95, 0.93, " 35.9 fb^{-1} (13 TeV)")

        text_latex = ROOT.TLatex()
        text_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignTop)
        text_latex.SetTextFont(42)
        text_latex.SetTextSize(0.03)
        start_y = 0.92
        diff_y = 0.07
        diff_y = 0.05
        for ind, line in enumerate(self.title.split('\n')):
            text_latex.DrawLatex(0.18, start_y - (ind*diff_y), line)

        # Do subplot
        if self.subplot_type and not self.subplot_type2:
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


            if self.xlim:
                self.subplot_container.GetXaxis().SetRangeUser(*self.xlim)

            if self.subplot_type == "ratio":
                self.subplot_container.SetMinimum(0)  # use this, not SetRangeUser()
                if self.subplot_limits:
                    self.subplot_container.SetMinimum(self.subplot_limits[0])  # use this, not SetRangeUser()
                    self.subplot_container.SetMaximum(self.subplot_limits[1])  # use this, not SetRangeUser()
                else:
                    # Make sure that the upper limit is the largest bin of the contributions,
                    # so long as it is within 1.5 and some upper limit
                    harrays = [cu.th1_to_arr(h) for h in self.subplot_contributions]
                    try:
                        bin_maxs = [np.max(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_maxs = [0]
                    except IndexError:
                        bin_maxs = [0]
                    self.subplot_container.SetMaximum(min(3, max(1.5, 1.2*max(bin_maxs))))

                    # Make sure the lower limit is the smallest bin of the contributions,
                    # so long as it is within 0 and 0.5
                    try:
                        bin_mins = [np.min(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_mins = [0]
                    except IndexError:
                        bin_mins = [1]
                    self.subplot_container.SetMinimum(min(0.5, 0.9*min(bin_mins)))

                xax = modifier.GetXaxis()
                self.subplot_line = ROOT.TLine(xax.GetXmin(), 1., xax.GetXmax(), 1.)
                if self.xlim:
                    self.subplot_line = ROOT.TLine(self.xlim[0], 1., self.xlim[1], 1.)
                self.subplot_line.SetLineStyle(2)
                self.subplot_line.SetLineWidth(2)
                self.subplot_line.SetLineColor(ROOT.kBlack)
                self.subplot_line.Draw()

            self._rescale_plot_labels(self.subplot_container, self.subplot_pad_height)
            xax = self.subplot_container.GetXaxis()
            xax.SetTitleOffset(self.subplot_container.GetXaxis().GetTitleOffset()*3)
            if xax.IsAlphanumeric():
                if xax.GetNbins() > 10:
                    xax.LabelsOption("v")
                    xax.SetLabelSize(0.08)
            self.subplot_container.GetYaxis().SetNdivisions(505)

            self.subplot_pad.Update()
            self.canvas.Update()

        elif self.subplot_type and self.subplot_type2:
            self._rescale_plot_labels(modifier, 1-self.subplot_pad_height-self.subplot_pad_height2)

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


            if self.xlim:
                self.subplot_container.GetXaxis().SetRangeUser(*self.xlim)

            if self.subplot_type == "ratio":
                self.subplot_container.SetMinimum(0)  # use this, not SetRangeUser()
                if self.subplot_limits:
                    self.subplot_container.SetMinimum(self.subplot_limits[0])  # use this, not SetRangeUser()
                    self.subplot_container.SetMaximum(self.subplot_limits[1])  # use this, not SetRangeUser()
                else:
                    # Make sure that the upper limit is the largest bin of the contributions,
                    # so long as it is within 1.5 and some upper limit
                    harrays = [cu.th1_to_arr(h) for h in self.subplot_contributions]
                    try:
                        bin_maxs = [np.max(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_maxs = [0]
                    except IndexError:
                        bin_maxs = [0]
                    self.subplot_container.SetMaximum(min(3, max(1.5, 1.2*max(bin_maxs))))

                    # Make sure the lower limit is the smallest bin of the contributions,
                    # so long as it is within 0 and 0.5
                    try:
                        bin_mins = [np.min(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_mins = [0]
                    except IndexError:
                        bin_mins = [1]
                    self.subplot_container.SetMinimum(min(0.5, 0.9*min(bin_mins)))

                xax = modifier.GetXaxis()
                self.subplot_line = ROOT.TLine(xax.GetXmin(), 1., xax.GetXmax(), 1.)
                if self.xlim:
                    self.subplot_line = ROOT.TLine(self.xlim[0], 1., self.xlim[1], 1.)
                self.subplot_line.SetLineStyle(2)
                self.subplot_line.SetLineWidth(2)
                self.subplot_line.SetLineColor(ROOT.kBlack)
                self.subplot_line.Draw()

            self._rescale_plot_labels(self.subplot_container, self.subplot_pad_height/0.7)
            xax = self.subplot_container.GetXaxis()
            xax.SetTitleOffset(9999999999)
            xax.SetLabelOffset(9999999999)
            # if xax.IsAlphanumeric():
            #     if xax.GetNbins() > 10:
            #         xax.LabelsOption("v")
            #         xax.SetLabelSize(0.08)
            self.subplot_container.GetYaxis().SetNdivisions(505)
            self.subplot_pad.Update()

            self.subplot_pad2.cd()
            self.subplot_container2.Draw(draw_opts)

            if self.subplot_title2 == None:
                if (self.subplot_type2 == "ratio"):
                    self.subplot_title2 = "#splitline{Ratio vs}{%s}" % (self.subplot2.label)
                elif (self.subplot_type2 == "diff"):
                    self.subplot_title2 = "#splitline{Difference}{vs %s}" % (self.subplot2.label)
                elif (self.subplot_type2 == "ddelta"):
                    self.subplot_title2 = "d#Delta/d#lambda"

            self.subplot_container2.SetTitle(";%s;%s" % (self.xtitle, self.subplot_title2))


            if self.xlim:
                self.subplot_container2.GetXaxis().SetRangeUser(*self.xlim)

            if self.subplot_type2 == "ratio":
                self.subplot_container2.SetMinimum(0)  # use this, not SetRangeUser()
                if self.subplot_limits2:
                    self.subplot_container2.SetMinimum(self.subplot_limits2[0])  # use this, not SetRangeUser()
                    self.subplot_container2.SetMaximum(self.subplot_limits2[1])  # use this, not SetRangeUser()
                else:
                    # Make sure that the upper limit is the largest bin of the contributions,
                    # so long as it is within 1.5 and some upper limit
                    harrays = [cu.th1_to_arr(h) for h in self.subplot_contributions2]
                    try:
                        bin_maxs = [np.max(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_maxs = [0]
                    except IndexError:
                        bin_maxs = [0]
                    self.subplot_container2.SetMaximum(min(3, max(1.5, 1.2*max(bin_maxs))))

                    # Make sure the lower limit is the smallest bin of the contributions,
                    # so long as it is within 0 and 0.5
                    try:
                        bin_mins = [np.min(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_mins = [0]
                    except IndexError:
                        bin_mins = [1]
                    self.subplot_container2.SetMinimum(min(0.5, 0.9*min(bin_mins)))

                xax = modifier.GetXaxis()
                self.subplot_line2 = ROOT.TLine(xax.GetXmin(), 1., xax.GetXmax(), 1.)
                if self.xlim:
                    self.subplot_line2 = ROOT.TLine(self.xlim[0], 1., self.xlim[1], 1.)
                self.subplot_line2.SetLineStyle(2)
                self.subplot_line2.SetLineWidth(2)
                self.subplot_line2.SetLineColor(ROOT.kBlack)
                self.subplot_line2.Draw()

            self._rescale_plot_labels(self.subplot_container2, self.subplot_pad_height2/0.8)
            xax = self.subplot_container2.GetXaxis()
            xax.SetTitleOffset(self.subplot_container2.GetXaxis().GetTitleOffset()*4)
            if xax.IsAlphanumeric():
                if xax.GetNbins() > 10:
                    xax.LabelsOption("v")
                    xax.SetLabelSize(0.08)
            self.subplot_container2.GetYaxis().SetNdivisions(505)

            self.subplot_pad2.Update()
            self.canvas.Update()


    def save(self, filename):
        """Save the plot to file. Do some check to make sure dir exists."""
        filename = os.path.abspath(filename)
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        self.canvas.SaveAs(filename)
        del self.canvas
