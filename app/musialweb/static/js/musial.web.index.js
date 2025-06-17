/**
 * Displays a notification element with a custom message to the user.
 *
 * @param {String} text The message to display.
 */
function displayNotification(text) {
  $("#menu").append(
    `<div class='notification'><i class="fa-duotone fa-spinner-third fa-spin fa-2xl"></i> ` +
      text +
      `</div>`
  );
}

/**
 * Removes the notification element from the DOM.
 */
function removeNotification() {
  $(".notification").remove();
}

/**
 * Displays a toast element with a custom error message to the user.
 *
 * @param {String} text The message to display.
 */
function displayError(text) {
  $("#menu").append(
    `<div class='alert'><i class="fa-duotone fa-circle-exclamation"></i> ` +
      text +
      `<button class='undecorated float-right' onclick='$(".alert").remove()'><i class="fa-solid fa-x"></i></button></div>`
  );
}

/**
 * Checks if the response status is 200.
 *
 * @param {Object} response Flask response object.
 * @returns {Boolean} True if the response status is 200, false otherwise.
 */
function responseOk(response) {
  return response.status === 200;
}

/**
 * Redirects the user to the about page.
 */
function redirectAbout() {
  window.open(window.location.origin + "/about", "_blank");
}

/**
 * Provides a blob as a downloadable file.
 *
 * @param {Blob} blob The file-like blob object whose content is stored.
 * @param {String} name The name of the file to store the specified blob.
 */
function downloadBlob(blob, name) {
  var download_link = document.createElement("a");
  download_link.href = window.URL.createObjectURL(new Blob([blob]));
  download_link.download = name;
  download_link.click();
  download_link.remove();
}
