class ProgressBar():
    """
    """
    TEMPLATE = "{0:<40s} | {1:>5.1f}% | {2:15s} [Evt rate : {3:>4.1f} kHz] [Eff : {4:8.4f}%]"
    def __init__(self, name):
        """

        """
        self.name = "[%s]" % name
        if len(self.name) >= 35:
            self.name = self.name[0:26] + "..." + self.name[-11:]
        self.bar = ""


    def update(self, counts, performance, phys_summary):
        self.bar = self.TEMPLATE.format(
                self.name,
                100. * (counts["completed"] / counts["all"]),
                self.get_bar(counts),
                0 if performance["time"] == 0 else 0.001 * (phys_summary["n_events_initial"] / performance["time"]),
                0 if phys_summary["n_events_initial"] == 0 else 100. * (phys_summary["n_events_selected"]["nominal"] / phys_summary["n_events_initial"])
        )
        

    def get_bar(self, counts):
        return "{0:>4d}/{1:<4d} done, {2:>4d} running.".format(counts["completed"], counts["all"], counts["running"])
