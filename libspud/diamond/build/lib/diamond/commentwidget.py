#    This file is part of Diamond.
#
#    Diamond is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Diamond is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Diamond.  If not, see <http://www.gnu.org/licenses/>.

import sys
from gi.repository import GObject as gobject
from gi.repository import Gtk
from gi.repository import Pango as pango

class CommentWidget(Gtk.Frame):
    __gsignals__ = {
        "on-store": (
            gobject.SignalFlags.RUN_LAST,
            gobject.TYPE_NONE,
            ()
        ),
    }
    fontsize = 12

    def __init__(self):
        super().__init__()

        # Create scrolled window and text view
        scrolled_window = Gtk.ScrolledWindow()
        scrolled_window.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)

        self.textView = Gtk.TextView()
        self.textView.set_editable(False)
        self.textView.set_wrap_mode(Gtk.WrapMode.WORD)
        self.textView.set_cursor_visible(False)
        self.textView.modify_font(pango.FontDescription(str(self.fontsize)))
        self.textView.connect("focus-in-event", self.focus_in)
        self.textView.connect("focus-out-event", self.focus_out)

        scrolled_window.add(self.textView)

        # Store buffer and create a persistent comment tag
        self.buffer = self.textView.get_buffer()
        self.comment_tag = self.buffer.create_tag(
            "comment_tag",
            foreground="grey"
        )

        # Frame label
        label = Gtk.Label()
        label.set_markup("<b>Comment</b>")

        self.set_shadow_type(Gtk.ShadowType.NONE)
        self.set_label_widget(label)
        self.add(scrolled_window)

        # State
        self.comment_tree = None
        self.interacted = False

    def update(self, node):
        """
        Update the widget with the given node
        """
        # Store previous content if needed
        self.store()

        # No active node
        if node is None or not node.active:
            self.buffer.set_text("")
            self.textView.set_cursor_visible(False)
            self.textView.set_editable(False)
            try:
                self.textView.set_tooltip_text("")
                self.textView.set_property("has-tooltip", False)
            except Exception:
                pass
            self.interacted = False
            return

        # Retrieve comment data
        self.comment_tree = comment_tree = node.get_comment()
        buf = self.buffer
        tag = self.comment_tag

        if comment_tree is None:
            buf.set_text("No comment")
            self.textView.set_cursor_visible(False)
            self.textView.set_editable(False)
            tag.set_property("foreground", "grey")
            try:
                self.textView.set_tooltip_text("")
                self.textView.set_property("has-tooltip", False)
            except Exception:
                pass
        else:
            text = comment_tree.data if comment_tree.data is not None else "(string)"
            buf.set_text(text)
            if node.active:
                self.textView.set_cursor_visible(True)
                self.textView.set_editable(True)
                tag.set_property("foreground", "black")
            else:
                self.textView.set_cursor_visible(False)
                self.textView.set_editable(False)
                tag.set_property("foreground", "grey")

        start, end = buf.get_bounds()
        buf.apply_tag(tag, start, end)
        self.interacted = False

    def store(self):
        """
        Store data in the node comment.
        """
        comment_tree = self.comment_tree
        if comment_tree is None or not self.interacted:
            return

        start, end = self.buffer.get_bounds()
        new_comment = self.buffer.get_text(start, end, True)

        if new_comment != comment_tree.data:
            if new_comment == "":
                comment_tree.data = None
                comment_tree.active = False
            else:
                comment_tree.set_data(new_comment)
                comment_tree.active = True
                self.emit("on-store")

    def focus_in(self, widget, event):
        """
        Called when the comment widget gains focus. Removes the placeholder.
        """
        comment_tree = self.comment_tree
        if comment_tree is not None and not self.interacted:
            self.interacted = True
            if comment_tree.data is None:
                self.buffer.set_text("")
        return False

    def focus_out(self, widget, event):
        """
        Called when the comment widget loses focus. Stores the comment.
        """
        self.store()
        return False

    def increase_font(self):
        self.fontsize += 2
        self.textView.modify_font(pango.FontDescription(str(self.fontsize)))

    def decrease_font(self):
        if self.fontsize > 0:
            self.fontsize -= 2
            self.textView.modify_font(pango.FontDescription(str(self.fontsize)))

# Register GObject type
gobject.type_register(CommentWidget)
