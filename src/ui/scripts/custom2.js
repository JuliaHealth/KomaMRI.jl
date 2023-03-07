if (window.module) module = window.module;

const { shell } = require('electron');

//Open links in new Window
mainWindow.webContents.on('new-window', function(e, url) {
    e.preventDefault();
    shell.openExternal(url);
});
