import AppKit

guard CommandLine.arguments.count == 3,
      let icon = NSImage(contentsOfFile: CommandLine.arguments[1]),
      NSWorkspace.shared.setIcon(icon, forFile: CommandLine.arguments[2], options: []) else {
    fatalError("PepSeMan launcher icon could not be set")
}
