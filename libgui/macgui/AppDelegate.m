//
//  AppDelegate.m
//  VagabondViewer
//
//  Created by Helen Ginn on 2017
//  Copyright (c) 2017 Helen Ginn. All rights reserved.
//

#import "AppDelegate.h"

@implementation AppDelegate

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
	int width = [[NSScreen mainScreen] frame].size.width;
	int height = [[NSScreen mainScreen] frame].size.height;

	[self.window setFrame:NSMakeRect(10, 10, width - 20, height - 20) display:YES];

	[self.window setCollectionBehavior:
	 NSWindowCollectionBehaviorFullScreenPrimary];

	[self.window makeFirstResponder:self.view];
	[self.window setAcceptsMouseMovedEvents:YES];
}

@end
