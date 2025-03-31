$(document).ready(function () {
    // Function to load announcements and insert them into the DOM
    function loadAnnouncements(url, targetSelector) {
        $.get(url, function (data) {
            // Create a new div and add the fetched content
            const newDiv = $("<div>").html(data);

            // Insert the new div after the target element
            $(targetSelector).after(newDiv);
        }).fail(function () {
            console.error("Failed to load the message from the external file: " + url);
        });
    }

    // Load Immcantation announcements and insert into a new div
    loadAnnouncements("https://raw.githubusercontent.com/immcantation/immcantation/refs/heads/master/announcements.html", "#welcome-to-the-immcantation-portal > h1");

});